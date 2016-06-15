package pmmh

import breeze.stats.distributions.Multinomial
import org.apache.log4j.{Level, Logger}

/**
  * Created by hshzcs on 6/7/16.
  */
object lv12_spark {
    import java.io.{Writer, File, PrintWriter}

    import breeze.linalg.DenseVector
    import breeze.stats.distributions.{Uniform, Exponential, Gaussian, Poisson}

    import org.apache.spark.SparkContext

    import spire.math.Fractional
    import spire.implicits._

    import scala.annotation.tailrec
    import scala.collection.mutable.ListBuffer
    import scala.io.Source
    import scala.math._
    import scala.math.log
    import scala.util.control.Breaks._

    type LvParameter2 = List[Double]
    type TS[O] = List[(Double, O)]

    def main(args: Array[String]): Unit = {
        val start = System.currentTimeMillis()
        println("Starting...")
        val its = if (args.length == 0) 1000 else args(0).toInt
        println("Running for " + its + " iters:")
        runModel(its)
        println("\nDone.")
        println("total time: " + (System.currentTimeMillis() - start) / 1000 + "s")
    }

    def runModel(its: Int) = {
        val rawData = Source.fromFile("/home/hshzcs/Desktop/scania-pmcmc/Spark/LV12data-4.txt").getLines.toList.take(6).map(_.split("\t").map(_.trim.toDouble))
        val data = rawData map { x => (x(0), List(x(1), x(2), x(3))) }

        val context = helper.Configuration.connectToSparkCluster("LV12")
        val rootLogger = Logger.getRootLogger()
        rootLogger.setLevel(Level.ERROR)

        val mll = pfPropPar_spark2(1000, simPrior2, 0.0, stepLV2, obsLik2, data, context)

        val s = new PrintWriter(new File("/home/hshzcs/Desktop/scania-pmcmc/Spark/mcmc-out-lv12.csv"))
        s.write("ll,th1,th2,th3,th4,th5,")
        s.write(((0 to 20 by 4) map { n => "prey" + n + ",pred1_" + n + ",pred2_" + n}).mkString(",") + "\n")

        val pmmhOutput = runPmmhPath2(s, its, List(10.0, 0.005, 0.0025, 6.0, 3.0), mll, peturb2)

        s.close
        context.stop()
    }

    def simPrior2(n: Int, t: Double, th: LvParameter2): Vector[Array[Int]] = {
        val prey = new Poisson(800.0).sample(n).toArray
        var predator = ListBuffer(prey)
        predator += new Poisson(500.0).sample(n).toArray
        predator += new Poisson(600.0).sample(n).toArray
        predator.toArray.transpose.toVector
    }

    def obsLik2(s: Array[Int], o: List[Double]): Double = {
        var result = 0.0
        for(i <- 0 to (o.size-1)){
            result += Gaussian(s(i), 100).logPdf(o(i))
        }
        result
    }

    @tailrec def stepLV2(xorig: Array[Int], t0: Double, dt: Double, th: List[Double]): Array[Int] = {
        val x = xorig.clone()
        if (dt <= 0.0) x
        else {
            val n = x.size - 1//number of predator
            val prey = x(0)
            val h = new Array[Double](1+2*n)//ListBuffer(th(0) * prey)
            h(0) = th(0) * prey
            for(i <- 1 to n)
                h(i) = (th(i) * prey * x(i))
            for(i <- 1 to n)
                h(i+n) = (th(i+n) * x(i))

            val h0 = adSum(h)
            val t = if (h0 < 1e-10 || x(0) > 1e6) 1e99 else new Exponential(h0).draw
            //println("t: " + t)
            if (t > dt) x
            else{
                val u = new Uniform(0, 1).draw
                var cumulate = h(0)
                if(u < (cumulate/h0)){
                    x(0) += 1
                    stepLV2(x, t0 + t, dt - t, th)
                }
                else{
                    var tmp = -1
                    breakable {
                        for(i <- 1 to (h.size-1)){
                            cumulate += h(i)
                            if(u < (cumulate/h0)){
                                tmp = i
                                break
                                System.out.println("sth wrong in the flow, this line should not be reached")
                            }
                        }
                    }
                    if(tmp <= n) {
                        x(0) -= 1
                        x(tmp) += 1
                        stepLV2(x, t0 + t, dt - t, th)
                    }
                    else{
                        x(tmp-n) -= 1
                        stepLV2(x, t0 + t, dt - t, th)
                    }
                }
            }
        }
    }

    def adSum(ad: Array[Double]) = {
        var sum = 0.0
        var i = 0
        while (i < ad.length) {
            sum += ad(i)
            i += 1
        }
        sum
    }

    def diff[T: Fractional](l: Iterable[T]): Iterable[T] = {
        (l.tail zip l) map { x => x._1 - x._2 }
    }

    def peturb2(th: List[Double]): List[Double] = {
        th map {x => x * exp(Gaussian(0, 0.014).draw)}
    }

    def sample(n: Int, prob: DenseVector[Double]): Vector[Int] = {
        Multinomial(prob).sample(n).toVector
    }

    def mean(it: Iterable[Double]): Double = {
        it.sum / it.size
    }

    def getCol2(a: Vector[List[Array[Int]]], n: Int): Vector[Array[Int]] = a map {_(n)}

    //Spark, multiple predators
    def pfPropPar_spark2(
                                n: Int,
                                simx0: (Int, Double, LvParameter2) => Vector[Array[Int]],
                                t0: Double,
                                stepFun: (Array[Int], Double, Double, LvParameter2) => Array[Int],
                                dataLik: (Array[Int], List[Double]) => Double,
                                data: TS[List[Double]],
                                context: SparkContext): (LvParameter2 => (Double, List[Array[Int]])) = {
        val (times, obs) = data.unzip
        //List(0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        //time difference between two prey-predator guesses
        val deltas = diff(t0 :: times)
        (th: LvParameter2) => {
            val x0 = simx0(n, t0, th)
            @tailrec def pf(ll: Double, x: Vector[List[Array[Int]]], t: Double, deltas: Iterable[Double], obs: List[List[Double]]): (Double, List[Array[Int]]) =
                if(obs.size > 0){
                    val head = obs.head
                    val xp = if (deltas.head == 0) getCol2(x, 0) else (context.parallelize(getCol2(x, 0), 240) map { l => stepFun(l, t, deltas.head, th) }).collect().toVector
                    val lw = xp map { l => dataLik(l, head) }
                    val max = lw.max
                    val w = lw map { x => exp(x - max) }
                    val rows = sample(n, DenseVector(w.toArray))
                    val xpp = rows map {r => xp(r) :: x(r)}
                    pf(ll + max + log(mean(w)), xpp, t + deltas.head, deltas.tail, obs.tail)
                }
                else{
                    (ll, x.head.reverse)
                }
            pf(0, x0 map { _ :: Nil }, t0, deltas, obs)
        }
    }

    def runPmmhPath2(s: Writer, iters: Int, initialState: List[Double], mll: List[Double] => (Double, List[Array[Int]]), peturb: List[Double] => List[Double]): List[List[Double]] = {
        @tailrec def pmmhAcc(itsLeft: Int, currentState: List[Double], currentMll: Double, currentPath: List[Array[Int]], allIts: List[List[Double]]): List[List[Double]] = {
            System.err.print(itsLeft.toString + " ")
            s.write(currentMll+",")
            s.write(currentState.mkString(",") + ",")
            s.write((currentPath.tail map { x => x.mkString(",")}).mkString(",") + "\n")//tail only applies to Spark version
            if (itsLeft == 0) allIts
            else {
                val prop = peturb(currentState)
                val (propMll, propPath) = mll(prop)//smc
                //s.write("propMll: " + propMll+"\n")
                //println("prop: " + prop)
                //println("propPath: " + (propPath.zip(0 to 30) map { x => x._2 + ": " + x._1.mkString(",")}).mkString("|"))
                //s.write("currentMll: " + currentMll+"\n")
                //println("currentParam: " + currentState)
                //println("currentPath: " + (currentPath.zip(0 to 30) map { x => x._2 + ": " + x._1.mkString(",")}).mkString("|"))
                //s.write("diff: " + (propMll - currentMll)+"\n")
                val tmp = log(Uniform(0, 1).draw)
                //s.write("likelihood: " + tmp+"\n")
                if (tmp < propMll - currentMll) {
                    //s.write("update accepted\n")
                    pmmhAcc(itsLeft - 1, prop, propMll, propPath, prop :: allIts)
                } else {
                    //s.write("update rejected\n")
                    pmmhAcc(itsLeft - 1, currentState, currentMll, currentPath, currentState :: allIts)
                }
            }
        }
        pmmhAcc(iters, initialState, -1e99, mll(initialState)._2, Nil).reverse
    }

}
