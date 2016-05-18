package pmmh

import scala.collection.mutable.ListBuffer

/**
  * Created by hshzcs on 5/18/16.
  */
object lvsim {

    import java.io.{File, PrintWriter}

    import breeze.stats.distributions._
    import pfilter._
    import pmmh._

    import scala.annotation.tailrec
    import scala.io.Source
    import scala.math.exp
    import scala.util.control.Breaks._

    //type Time = Double

    //type LogLik = Double

    type LvParameter2 = List[Double]

    type TS[O] = List[(Double, O)]

    def simPrior(n: Int, t: Double, th: LvParameter): Vector[(Int, Int)] = {
        val prey = new Poisson(50.0).sample(n).toVector
        val predator = new Poisson(100.0).sample(n).toVector
        prey.zip(predator) map { x => (x._1, x._2) }
    }

    def simPrior2(n: Int, t: Double, th: LvParameter2): Vector[Array[Int]] = {
        val prey = new Poisson(800.0).sample(n).toArray
        var predator = ListBuffer(prey)
        for(i <- 1 to (th.size-1)/2)
            predator += new Uniform(0.0, 1000.0).sample(n).toArray map {_.toInt}
        predator.toArray.transpose.toVector
    }

    def obsLik(s: (Int, Int), o: Double): Double = {
        Gaussian(s._1, 10.0).logPdf(o)
    }

    def obsLik(s: (Int, Int), o: Array[Double]): Double = {
        Gaussian(s._1, 10.0).logPdf(o(0)) + Gaussian(s._2, 10.0).logPdf(o(1))
    }

    def obsLik2(s: Array[Int], o: Double): Double = {
        Gaussian(s(0), 10.0).logPdf(o)
    }


    //x0 = List(50, 100, 80, 40)

    //th = List(1.0, 0.005, 0.0025, 0.006, 0.6, 0.3, 0.72)
    def runModel(its: Int) = {
        val rawData = Source.fromFile("/home/hshzcs/Desktop/scania-pmcmc/Spark/LV13data2.txt").getLines.toList.take(31).map(_.split("\t").map(_.trim.toDouble))
        val data = rawData map { x => (x(0), x(1)) }//(rawData.indices.filter(_ % 2 == 0).map(rawData(_))).toList map { x => (x(0), x(1)) }
        //val rawData = Source.fromFile("/home/hshzcs/thesiscode/SparkApp/LVpreyNoise10.txt").getLines//.map(_.split(",").map(_.trim.toDouble))
        //val data2 = ((0 to 30 by 2).toList zip rawData.toList) map { x => (x._1.toDouble, x._2) }
        //parallel-spark
        //val context = helpers.Configuration.connectToSparkCluster("LV-v1")
        //val mll = pfPropPar_spark(100000, simPrior, 0.0, stepLV, obsLik, data, context)
        //parallel-spark2
        val context = helper.Configuration.connectToSparkCluster("LV13")
        val mll = pfPropPar_spark2(1000, simPrior2, 0.0, stepLV2, obsLik2, data, context)
        //parallel-scala
        //val mll = pfPropPar_scala(1000, simPrior, 0.0, stepLV, obsLik, data)
        //parallel-scala2
        //val mll = pfPropPar_scala2(1000, simPrior2, 0.0, stepLV2, obsLik2, data)
        //sequential
        //val mll = pfProp(100, simPrior, 0.0, stepLV, obsLik, data)
        val s = new PrintWriter(new File("/home/hshzcs/Desktop/scania-pmcmc/Spark/mcmc-out.csv"))
        // val s=new OutputStreamWriter(System.out)
        s.write("th1,th2,th3,th4,th5,th6,th7,")
        s.write(((0 to 30) map { n => "prey" + n + ",pred1_" + n + ",pred2_" + n +",pred3_" + n}).mkString(",") + "\n")
        //val pmmhOutput = runPmmhPath(s, its, new LvParameter(1.0, 0.005, 0.6), mll, peturb)
        val pmmhOutput = runPmmhPath2(s, its, List(10.0, 0.005, 0.0025, 0.006, 6.0, 3.0, 7.2), mll, peturb2)
        s.close
        //context.stop()
    }

    //th0: alpha, the reproduce rate of the prey
    //th1: beta,  the reproduce rate of the predator by consuming the prey
    //th2: gamma, the die rate of the predator
    //t0: current time
    //dt: time difference to next state
    //Gillespie algorithm
    @tailrec def stepLV(x: (Int, Int), t0: Double, dt: Double, th: LvParameter): (Int, Int) = {
        if (dt <= 0.0) x
        else {
            val h = (th.th0 * x._1, th.th1 * x._2 * x._1, th.th2 * x._2)
            val h0 = h._1 + h._2 + h._3
            val t = if (h0 < 1e-10 || x._1 > 1e6) 1e99 else new Exponential(h0).draw
            if (t > dt) x
            else {
                //val u = Random.nextDouble()
                val u = new Uniform(0, 1).draw // use simpler function!
                if (u < h._1 / h0)
                    stepLV((x._1 + 1, x._2), t0 + t, dt - t, th)
                else {
                    if (u < (h._1 + h._2) / h0)
                        stepLV((x._1 - 1, x._2 + 1), t0 + t, dt - t, th)
                    else
                        stepLV((x._1, x._2 - 1), t0 + t, dt - t, th)
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

    @tailrec def stepLV3(x: (Int, Int), t0: Double, dt: Double, th: LvParameter): (Int, Int) = {
        //println("dt"+dt)
        if (dt <= 0.0) x
        else {
            val h = (th.th0 * x._1, th.th1 * x._2 * x._1, th.th2 * x._2)
            val h0 = h._1 + h._2 + h._3
            val t = if (h0 < 1e-10 || x._1 > 1e6) 1e99 else new Exponential(h0).draw
            //println("t: "+t)
            if (t > dt) x
            else {
                //val u = Random.nextDouble()
                val u = new Uniform(0, 1).draw // use simpler function!
                if (u < h._1 / h0){
                    println((t0+t)+","+(x._1 + 1)+","+(x._2))
                    stepLV3((x._1 + 1, x._2), t0 + t, dt - t, th)
                }
                else {
                    if (u < (h._1 + h._2) / h0){
                        println((t0+t)+","+(x._1 - 1)+","+(x._2+1))
                        stepLV3((x._1 - 1, x._2 + 1), t0 + t, dt - t, th)
                    }
                    else{
                        println((t0+t)+","+(x._1)+","+(x._2 - 1))
                        stepLV3((x._1, x._2 - 1), t0 + t, dt - t, th)
                    }
                }
            }
        }
    }

    def peturb(th: LvParameter): LvParameter = {
        new LvParameter(th.th0 * exp(Gaussian(0, 0.01).draw), th.th1 * exp(Gaussian(0, 0.01).draw), th.th2 * exp(Gaussian(0, 0.01).draw))
    }

    def peturb2(th: List[Double]): List[Double] = {
        th map {x => x * exp(Gaussian(0, 0.01).draw)}
    }
}
