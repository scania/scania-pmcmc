package pmmh

/**
  * Created by hshzcs on 5/18/16.
  */
object pmmh {

    import java.io.Writer

    import breeze.stats.distributions._

    import scala.annotation.tailrec
    import scala.math.log

    //lv11
    def runPmmhPath(s: Writer, iters: Int, initialState: LvParameter, mll: LvParameter => (Double, List[(Int, Int)]), peturb: LvParameter => LvParameter): List[LvParameter] = {
        @tailrec def pmmhAcc(itsLeft: Int, currentState: LvParameter, currentMll: Double, currentPath: List[(Int, Int)], allIts: List[LvParameter]): List[LvParameter] = {
            System.err.print(itsLeft.toString + " ")
            s.write(currentState.toString + ",")
            s.write((currentPath map { x => x._1 + "," + x._2 }).mkString(",") + "\n")
            if (itsLeft == 0) allIts
            else {
                val prop = peturb(currentState)
                val (propMll, propPath) = mll(prop)//smc
                println("propMll: " + propMll+"\n")
                println("currentMll: " + currentMll+"\n")
                println("diff: " + (propMll - currentMll))
                val tmp = log(Uniform(0, 1).draw)
                println("likelihood: " + tmp)
                if (tmp < propMll - currentMll) {
                    println("update accepted\n")
                    pmmhAcc(itsLeft - 1, prop, propMll, propPath, prop :: allIts)
                } else {
                    println("update rejected\n")
                    pmmhAcc(itsLeft - 1, currentState, currentMll, currentPath, currentState :: allIts)
                }
            }
        }
        pmmhAcc(iters, initialState, (-1e99).toDouble, mll(initialState)._2, Nil).reverse
    }

    //lv13
    def runPmmhPath2(s: Writer, iters: Int, initialState: List[Double], mll: List[Double] => (Double, List[Array[Int]]), peturb: List[Double] => List[Double]): List[List[Double]] = {
        @tailrec def pmmhAcc(itsLeft: Int, currentState: List[Double], currentMll: Double, currentPath: List[Array[Int]], allIts: List[List[Double]]): List[List[Double]] = {
            System.err.print(itsLeft.toString + " ")
            s.write(currentState.mkString(",") + ",")
            s.write((currentPath map { x => x.mkString(",")}).mkString(",") + "\n")//tail only applies to Spark version
            if (itsLeft == 0) allIts
            else {
                val prop = peturb(currentState)
                val (propMll, propPath) = mll(prop)//smc
                s.write("propMll: " + propMll+"\n")
                //println("prop: " + prop)
                //println("propPath: " + (propPath.zip(0 to 30) map { x => x._2 + ": " + x._1.mkString(",")}).mkString("|"))
                s.write("currentMll: " + currentMll+"\n")
                //println("currentParam: " + currentState)
                //println("currentPath: " + (currentPath.zip(0 to 30) map { x => x._2 + ": " + x._1.mkString(",")}).mkString("|"))
                s.write("diff: " + (propMll - currentMll)+"\n")
                val tmp = log(Uniform(0, 1).draw)
                s.write("likelihood: " + tmp+"\n")
                if (tmp < propMll - currentMll) {
                    s.write("update accepted\n")
                    pmmhAcc(itsLeft - 1, prop, propMll, propPath, prop :: allIts)
                } else {
                    s.write("update rejected\n")
                    pmmhAcc(itsLeft - 1, currentState, currentMll, currentPath, currentState :: allIts)
                }
            }
        }
        pmmhAcc(iters, initialState, (-1e99).toDouble, mll(initialState)._2, Nil).reverse
    }
}
