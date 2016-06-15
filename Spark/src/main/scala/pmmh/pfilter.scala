package pmmh

import org.apache.spark.SparkContext
import scala.collection.parallel.immutable.{ParVector, ParSeq}

/**
  * Created by hshzcs on 5/18/16.
  */
object pfilter {

    import breeze.linalg.DenseVector
    import breeze.stats.distributions.Multinomial
    import lvsim._
    import spire.implicits._
    import spire.math._

    import scala.annotation.tailrec

    def diff[T: Fractional](l: Iterable[T]): Iterable[T] = {
        (l.tail zip l) map { x => x._1 - x._2 }
    }

    def sample(n: Int, prob: DenseVector[Double]): Vector[Int] = {
        Multinomial(prob).sample(n).toVector
    }

    def mean(it: Iterable[Double]): Double = {
        it.sum / it.size
    }

    //sequential, single predator
    def pfProp(
                      n: Int,
                      simx0: (Int, Double, LvParameter) => Vector[(Int, Int)],
                      t0: Double,
                      stepFun: ((Int, Int), Double, Double, LvParameter) => (Int, Int),
                      dataLik: ((Int, Int), Array[Double]) => Double,
                      data: TS[Array[Double]]): (LvParameter => (Double, List[(Int, Int)])) = {
        val (times, obs) = data.unzip
        val deltas = diff(t0 :: times)
        (th: LvParameter) => {
            val x0 = simx0(n, t0, th) // original, randomly generated n pair of prey and predator
            @tailrec def pf(ll: Double, x: Vector[List[(Int, Int)]], t: Double, deltas: Iterable[Double], obs: List[Array[Double]]): (Double, List[(Int, Int)]) =
                obs match {
                    case Nil => (ll, x(0).reverse)
                    case head :: tail => {
                        //println(x(0).size)
                        val xp = if (deltas.head == 0) x else (x map { l => stepFun(l.head, t, deltas.head, th) :: l })
                        //println("x: " + x)
                        val lw = xp map { l => dataLik(l.head, head) }
                        //println("lw: " + lw)
                        val max = lw.max
                        //println("max: " + max)
                        val w = lw map { x => exp(x - max) }
                        //println("w: " + w)
                        val rows = sample(n, DenseVector(w.toArray))
                        //println("rows: " + rows)
                        val xpp = rows map {
                            xp(_)
                        }
                        //println("xpp: " + xpp)
                        pf(ll + max + log(mean(w)), xpp, t + deltas.head, deltas.tail, tail)
                    }
                }
            pf(0, x0 map { _ :: Nil }, t0, deltas, obs)
        }
    }

    //sequential, multiple predators
    def pfProp2(
                                n: Int,
                                simx0: (Int, Double, LvParameter2) => Vector[Array[Int]],
                                t0: Double,
                                stepFun: (Array[Int], Double, Double, LvParameter2) => Array[Int],
                                dataLik: (Array[Int], List[Double]) => Double,
                                data: TS[List[Double]]): (LvParameter2 => (Double, List[Array[Int]])) = {
        val (times, obs) = data.unzip
        //List(0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        //time difference between two prey-predator guesses
        val deltas = diff(t0 :: times)
        (th: LvParameter2) => {
            val x0 = simx0(n, t0, th)
            @tailrec def pf(ll: Double, x: Vector[List[Array[Int]]], t: Double, deltas: Iterable[Double], obs: List[List[Double]]): (Double, List[Array[Int]]) =
            //obs match {
            //case Nil => (ll, x(0).reverse)
            //case head :: tail => {
                if(obs.size > 0) {
                    //println("ll: "+ ll)
                    val head = obs.head
                    val xp = if (deltas.head == 0) x else (x map { l => stepFun(l.head, t, deltas.head, th) :: l })
                    val lw = xp map { l => dataLik(l.head, head) }
                    val max = lw.max
                    val w = lw map { x => exp(x - max) }
                    val rows = sample(n, DenseVector(w.toArray))
                    val xpp = rows map { xp(_) } // the _th element of xp
                    pf(ll + max + log(mean(w)), xpp, t + deltas.head, deltas.tail, obs.tail)
                }
                else{
                    //println("ll: "+ ll)
                    (ll, x.head.reverse)
                }
            pf(0, x0 map { _ :: Nil }, t0, deltas, obs)
        }
    }

    def mean(it: ParSeq[Double]): Double = {
        it.sum / it.size
    }

    //parallel; single predator
    def pfPropPar_scala(
                         n: Int,
                         simx0: (Int, Double, LvParameter) => Vector[(Int, Int)],
                         t0: Double,
                         stepFun: ((Int, Int), Double, Double, LvParameter) => (Int, Int),
                         dataLik: ((Int, Int), Double) => Double,
                         data: TS[Double]): (LvParameter => (Double, List[(Int, Int)])) = {
        val (times, obs) = data.unzip
        //List(0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        //time difference between two prey-predator guesses
        val deltas = diff(t0 :: times)
        (th: LvParameter) => {
            val x0 = simx0(n, t0, th).par
            @tailrec def pf(ll: Double, x: ParVector[List[(Int, Int)]], t: Double, deltas: Iterable[Double], obs: List[Double]): (Double, List[(Int, Int)]) =
                obs match {
                    case Nil => {
                        println("ll: "+ll)
                        (ll, x(0).reverse)
                    }
                    case head :: tail => {
                        println("ll: "+ll)
                        val xp = if (deltas.head == 0) x else (x map { l => stepFun(l.head, t, deltas.head, th) :: l })
                        val lw = xp map { l => dataLik(l.head, head) }
                        val max = lw.max
                        val w = lw map { x => exp(x - max) }
                        val rows = sample(n, DenseVector(w.toArray)).par
                        val xpp = rows map { xp(_) } // the _th element of xp
                        pf(ll + max + log(mean(w)), xpp, t + deltas.head, deltas.tail, tail)
                    }
                }
            pf(0, x0 map { _ :: Nil }, t0, deltas, obs)
        }
    }

    //parallel; multiple predators
    def pfPropPar_scala2(
                               n: Int,
                               simx0: (Int, Double, LvParameter2) => Vector[Array[Int]],
                               t0: Double,
                               stepFun: (Array[Int], Double, Double, LvParameter2) => Array[Int],
                               dataLik: (Array[Int], List[Double]) => Double,
                               data: TS[List[Double]]): (LvParameter2 => (Double, List[Array[Int]])) = {
        val (times, obs) = data.unzip
        //List(0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        //time difference between two prey-predator guesses
        val deltas = diff(t0 :: times)
        (th: LvParameter2) => {
            val x0 = simx0(n, t0, th).par
            @tailrec def pf(ll: Double, x: ParVector[List[Array[Int]]], t: Double, deltas: Iterable[Double], obs: List[List[Double]]): (Double, List[Array[Int]]) =
                //obs match {
                    //case Nil => (ll, x(0).reverse)
                    //case head :: tail => {
                if(obs.size > 0) {
                    //println("ll: "+ ll)
                    val head = obs.head
                    val xp = if (deltas.head == 0) x else (x map { l => stepFun(l.head, t, deltas.head, th) :: l })
                    val lw = xp map { l => dataLik(l.head, head) }
                    val max = lw.max
                    val w = lw map { x => exp(x - max) }
                    val rows = sample(n, DenseVector(w.toArray)).par
                    val xpp = rows map { xp(_) } // the _th element of xp
                    pf(ll + max + log(mean(w)), xpp, t + deltas.head, deltas.tail, obs.tail)
                }
                else{
                    //println("ll: "+ ll)
                    (ll, x.head.reverse)
                }
            pf(0, x0 map { _ :: Nil }, t0, deltas, obs)
        }
    }

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
                //println("ll: "+ ll)
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
                //println("ll: "+ ll)
                (ll, x.head.reverse)
            }
            pf(0, x0 map { _ :: Nil }, t0, deltas, obs)
        }
    }

    //Spark; single predator
    def pfPropPar_spark(
                         n: Int,
                         simx0: (Int, Double, LvParameter) => Vector[(Int, Int)],
                         t0: Double,
                         stepFun: ((Int, Int), Double, Double, LvParameter) => (Int, Int),
                         dataLik: ((Int, Int), Double) => Double,
                         data: TS[Double],
                         context: SparkContext): (LvParameter => (Double, List[(Int, Int)])) = {
        val (times, obs) = data.unzip
        val deltas = diff(t0 :: times)
        (th: LvParameter) => {
            val x0 = simx0(n, t0, th) //.par
            @tailrec def pf(ll: Double, x: Vector[List[(Int, Int)]], t: Double, deltas: Iterable[Double], obs: List[Double]): (Double, List[(Int, Int)]) =
                if(obs.size > 0){
                        //val subx = context.parallelize(x.grouped(200).toVector, 16)//val subx = x.grouped(10).toVector
                        val head = obs.head
                        val xp = if (deltas.head == 0) getCol(x, 0) else (context.parallelize(getCol(x, 0), 12) map { l => stepFun(l, t, deltas.head, th) }).collect().toVector
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
            println("n: " + n)
            pf(0, x0 map { _ :: Nil}, t0, deltas, obs)
        }
    }

    def getCol(a: Vector[List[(Int, Int)]], n: Int): Vector[(Int, Int)] = a map {_(n)}

    def getCol2(a: Vector[List[Array[Int]]], n: Int): Vector[Array[Int]] = a map {_(n)}
}
