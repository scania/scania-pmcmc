package pmmh

/**
  * Created by hshzcs on 5/18/16.
  */
//@SerialVersionUID(-4021364619531676066L)
class LvParameter(th0x: Double, th1x: Double, th2x: Double) extends Serializable {
    val th0 = th0x
    val th1 = th1x
    val th2 = th2x

    override def toString = "" + th0 + "," + th1 + "," + th2
}
