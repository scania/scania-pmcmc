package pmmh

import breeze.stats.distributions.Gaussian

import scala.math._

/**
  * Created by hshzcs on 5/18/16.
  */
object bayeskit {

    def main(args: Array[String]): Unit = {
        val start = System.currentTimeMillis()
        println("Starting...")
        val its = if (args.length == 0) 100 else args(0).toInt
        println("Running for " + its + " iters:")
        lvsim.runModel(its)
        println("\nDone.")
        println("total time: " + (System.currentTimeMillis() - start)/1000 + "s")
        //var x0 = List(100, 150, 80)

        //(1.0, 0.005, 0.6)

        //val ths = simPrior2(1000, List(5.2, 0.0147, 0.014, 0.0125, 3.1, 2.95, 2.6))
        //val th = List(1.0, 0.005, 0.002, 0.6, 0.4)

        //for(i <- 0 to ths.size - 1){
            //println(ths(i))
        //    for(j <- (0 to 29) ) {
        //        print("x" + (j+1) + ": ")
        //      x0 = lvsim.stepLV2(x0, j, j+1, th)
        //      println(x0)
        //  }
            //println(ths(i))
        //if(x0(0) != 0 && x0(0) != 1000001) {
        //      println(ths(i))
        //      println(x0)
        //      //println("!!!!!!!!!!!!!!!")
        //  }
            //x0 = List(150, 60, 80, 70)
        //}
        //for(i <- 0 to 0)
        //getValidList()
        //println(System.getProperty("java.library.path"))
    }

    //generated lv13data
    def getValidList() = {
        var x0 = Array(800,500,600)

        var tmp = Array(0, 0, 0)

        val th = List(10.0, 0.005, 0.0025, 6.0, 3.0)

        //var result = ListBuffer(x0)

        var j=0

        //println("x0: "+x0.mkString(","))

        while (j < 1000){
            //print("x" + (j+1) + ": ")
            tmp = lvsim.stepLV2(x0, 0, 4, th)
            //x0 = lvsim.stepLV3(x0, j, 1, new LvParameter(1.0, 0.005, 0.6))
            println(tmp.mkString(","))
            //println((j+1)+","+x0._1+","+x0._2)
            j = j+1
        }
    }

    def getValidList_lv11() = {
        var x0 = Array(35, 100)

        var tmp = Array(0, 0)

        val th = List(1.0, 0.005, 0.6)

        //var result = ListBuffer(x0)

        var j=0

        println(x0.mkString(","))

        while (j < 100){
            //print("x" + (j+1) + ": ")
            x0 = lvsim.stepLV2(x0, 0, 1, th)
            //x0 = lvsim.stepLV3(x0, j, 1, new LvParameter(1.0, 0.005, 0.6))
            println(x0.mkString(","))
            //println((j+1)+","+x0._1+","+x0._2)
            j = j+1
        }
    }
}
