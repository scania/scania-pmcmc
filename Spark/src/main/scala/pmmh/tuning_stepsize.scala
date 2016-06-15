package pmmh

import java.io.{Writer, File, PrintWriter}

import breeze.numerics.exp
import breeze.stats.distributions.{Uniform, Gaussian}
import org.apache.spark.SparkContext
import pfilter._
import lvsim._

import scala.annotation.tailrec
import scala.io.Source
import scala.math._

/**
  * Created by hshzcs on 5/31/16.
  */
/*
* result:
*
* sigma     acceptance rate
* 0.02      0.15    0.1     0.13    0.15    0.14    0.13    0.16    0.15    0.16    0.11
* 0.015     0.27    0.19    0.15    0.22    0.18    0.19    0.14    0.19    0.22    0.15    0.26    0.22    0.17    0.15    0.19    0.17    0.24    0.18    0.17    0.19
* 0.014     0.23    0.29    0.26    0.16    0.28    0.22    0.35    0.16    0.27    0.28    0.18    0.25    0.28    0.22    0.23    0.2     0.18    0.21    0.19    0.23
* 0.013     0.24    0.23    0.27    0.31    0.23    0.3     0.2     0.3     0.24    0.19
* 0.012     0.25    0.22    0.24    0.29    0.22    0.32    0.31    0.24    0.31    0.33
* 0.01      0.25    0.29    0.27    0.4     0.26    0.33    0.36    0.25    0.29    0.32
*
* average:
* 0.02:     0.138
* 0.015:    0.192
* 0.014:    0.2335
* 0.013:    0.251
* 0.012:    0.273
* 0.01:     0.302
*
* */
object tuning_stepsize {

    def main(args: Array[String]): Unit = {
        val start = System.currentTimeMillis()
        println("Starting...")
        //val sigma = 0.05
        val seq = List(0.014, 0.014, 0.014, 0.014, 0.014)
        //val context = helper.Configuration.connectToSparkCluster("LV12-tuning-sigma="+sigma)
        //val ar = getAcceptanceRate(100, sigma, context)
        //println("sigma: " + sigma + ", acceptance rate: ")
        val context = helper.Configuration.connectToSparkCluster("LV12-tuning")
        val s = new PrintWriter(new File("/home/hshzcs/Desktop/scania-pmcmc/Spark/tuning.csv"))
        for (x <- seq) {
            println("sigma: " + x)
            s.write("sigma: " + x+"\n")
            val ar = getAcceptanceRate(100, x, context)
            println("acceptance rate: " + ar + "\n")
            s.write("acceptance rate: " + ar +"\n\n")
        }
        s.close()
        println("\nDone.")
        println("total time: " + (System.currentTimeMillis() - start) / 1000 + "s")
    }

    def getAcceptanceRate(its: Int, sigma: Double, context: SparkContext) = {
        val rawData = Source.fromFile("/home/hshzcs/Desktop/scania-pmcmc/Spark/LV12data-4.txt").getLines.toList.take(6).map(_.split("\t").map(_.trim.toDouble))
        val data = rawData map { x => (x(0),List(x(1), x(2), x(3))) }//(rawData.indices.filter(_ % 2 == 0).map(rawData(_))).toList map { x => (x(0), x(1)) }
        //val rawData = Source.fromFile("/home/hshzcs/Desktop/scania-pmcmc/Spark/LVpreyNoise10.txt").getLines//.map(_.split(",").map(_.trim.toDouble))
        //val data = ((0 to 30 by 2).toList zip rawData.toList) map { x => (x._1.toDouble, x._2.toDouble) }
        //parallel-spark
        //val context = helper.Configuration.connectToSparkCluster("LV11")
        //val mll = pfPropPar_spark(500000, simPrior, 0.0, stepLV, obsLik, data, context)
        //parallel-spark2
        //val context = helper.Configuration.connectToSparkCluster("LV12-tuning-sigma="+sigma)
        val mll = pfPropPar_spark2(1000, simPrior2, 0.0, stepLV2, lv12_spark.obsLik2, data, context)
        //parallel-scala
        //val mll = pfPropPar_scala(1000, simPrior, 0.0, stepLV, obsLik, data)
        //parallel-scala2
        //val mll = pfPropPar_scala2(1000, simPrior2, 0.0, stepLV2, obsLik2, data)
        //sequential
        //val mll = pfProp(100, simPrior, 0.0, stepLV, obsLik, data)
        val s = new PrintWriter(new File("/home/hshzcs/Desktop/scania-pmcmc/Spark/mcmc-out-lv12-tuning-"+sigma+".csv"))
        // val s=new OutputStreamWriter(System.out)
        //s.write("th1,th2,th3,")
        //s.write(((0 to 30 by 2) map { n => "prey" + n + ",pred1_" + n + ",pred2_" + n}).mkString(",") + "\n")
        s.write("ll,th1,th2,th3,th4,th5,")
        s.write(((0 to 20 by 4) map { n => "prey" + n + ",pred1_" + n + ",pred2_" + n}).mkString(",") + "\n")
        //val pmmhOutput = runPmmhPath(s, its, new LvParameter(1.0, 0.005, 0.6), mll, peturb)
        val pmmhOutput = runPmmhPath(s, its, List(10.0, 0.005, 0.0025, 6.0, 3.0), mll, sigma)
        //val pmmhOutput = runPmmhPath2(s, its, List(1.0, 0.005, 0.0025, 0.6, 0.3), mll, peturb2)
        s.close
        //context.stop()
        pmmhOutput
    }

    def peturb(th: List[Double], sigma: Double): List[Double] = {
        th map {x => x * exp(Gaussian(0, sigma).draw)}
    }

    //lv13
    def runPmmhPath(s: Writer, iters: Int, initialState: List[Double], mll: List[Double] => (Double, List[Array[Int]]), sigma: Double): Double = {
        @tailrec def pmmhAcc(itsLeft: Int, currentState: List[Double], currentMll: Double, currentPath: List[Array[Int]], acceptedNum: Int): Double = {
            System.err.print(itsLeft.toString + " ")
            s.write(currentMll+",")
            s.write(currentState.mkString(",") + ",")
            s.write((currentPath.tail map { x => x.mkString(",")}).mkString(",") + "\n")//tail only applies to Spark version
            if (itsLeft == 0) acceptedNum/(iters+0.0)
            else {
                val prop = peturb(currentState, sigma)
                val (propMll, propPath) = mll(prop)//smc

                s.write("propMll: " + propMll+"\n")
                s.write("currentMll: " + currentMll+"\n")
                s.write("diff: " + (propMll - currentMll)+"\n")

                val tmp = log(Uniform(0, 1).draw)
                if (tmp < propMll - currentMll) {
                    s.write("update accepted\n")
                    pmmhAcc(itsLeft - 1, prop, propMll, propPath, acceptedNum + 1)
                } else {
                    s.write("update rejected\n")
                    pmmhAcc(itsLeft - 1, currentState, currentMll, currentPath, acceptedNum)
                }
            }
        }
        pmmhAcc(iters, initialState, (-1e99).toDouble, mll(initialState)._2, 0)
    }
}
