package helper

import java.io.File
import java.util.Map.Entry

import com.typesafe.config.ConfigValueType._
import com.typesafe.config.{ConfigFactory, ConfigValue}
import org.apache.spark.{SparkConf, SparkContext}

import scala.util.Try

object Configuration {

    // Initialize most of the variables

    val confPath = System.getProperty("os.name") match {
        case "Mac OS X" => "/path/to/conf/"
        case "Linux" => "/home/hshzcs/Desktop/scania-pmcmc/Spark/"
    }

    val sparkConfPath = "/cluster/sesonas13/hshzcs/spark"
    val sparkHome = "/share/apps/spark/spark-1.6.0"
    //val sparkMyHome = "/home/hshzcs/Applications/spark-1.5.2/"


    private def getPropertiesList(key: String) = {
        val list = Try(config.getConfig(key).entrySet().toArray).getOrElse(Array())
        list.map(x => {
            val p = x.asInstanceOf[Entry[String, ConfigValue]]
            val k = p.getKey

            val v = p.getValue.valueType match {
                case BOOLEAN => config.getBoolean(key + "." + k)
                case STRING => config.getString(key + "." + k)
                case NUMBER => config.getDouble(key + "." + k)
                case _ => config.getString(key + "." + k)
            }
            (k.replace("_", "."), v.toString)
        })
    }

    // Spark configuration - loading from "resources/application.conf"
    private val config = ConfigFactory.load()
    lazy val SPARK_MASTER_HOST = Try(config.getString("spark.master_host")).getOrElse("local")
    lazy val SPARK_MASTER_PORT = Try(config.getInt("spark.master_port")).getOrElse(7077)
    lazy val SPARK_HOME = Try(config.getString("spark.home")).getOrElse(sparkConfPath)
    lazy val SPARK_MEMORY = Try(config.getString("spark.memory")).getOrElse("1g")
    lazy val SPARK_DRIVER_MEMORY = Try(config.getString("spark.driver.memory")).getOrElse("1g")
    lazy val SPARK_OPTIONS = getPropertiesList("spark.options")
    lazy val SPARK_DEFAULT_PAR = Try(config.getString("spark_default_parallelism")).getOrElse("8")

    def connectToSparkCluster(appName: String): SparkContext = {
        // get the name of the packaged
        val thisPackagedJar = new File("/home/hshzcs/Desktop/scania-pmcmc/Spark/target/scala-2.10").listFiles.
                filter(x => x.isFile && x.getName.toLowerCase.takeRight
                (4) == ".jar").toList.map(_.toString)

        // Scan for external libraries in folder 'lib'
        // All the JAR files in this folder will be shipped to the cluster
        /*val libs = new File("/home/hshzcs/Desktop/scania-pmcmc/Spark/lib").listFiles.
                filter(x => x.isFile && x.getName.toLowerCase.takeRight(4) == ".jar").
                toList.map(_.toString)*/

        val master = "spark://10.6.3.238:7077"

        // Spark Context configuration
        val scConf =
            SPARK_OPTIONS.fold(
                new SparkConf()
                        .setMaster(master)
                        .setAppName(appName)
                        .set("spark.executor.memory", SPARK_MEMORY)
                        .setSparkHome(sparkHome)
                        .setJars(thisPackagedJar)
                        .set("spark.default.parallelism", SPARK_DEFAULT_PAR)
                        .set("spark.eventLog.enabled", "true")
                        .set("spark.eventLog.dir", sparkConfPath + "/workerLog")
                        //          .set("spark.task.cpus","8")
                        .set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
            )((c, p) => {
                // apply each spark option from the configuration file in section "spark.options"
                val (k, v) = p.asInstanceOf[(String, String)]
                c.asInstanceOf[SparkConf].set(k, v)
            }
            ).asInstanceOf[SparkConf]
        // Create and return the spark context to be used through the entire KYC application
        new SparkContext(scConf)
    }
}
