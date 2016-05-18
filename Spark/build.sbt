import scala.util.Properties

lazy val versions = Map(
    'slf4j -> Properties.envOrElse("SL4J_VERSION", "1.7.12"),
    'spark -> Properties.envOrElse("SPARK_VERSION", "1.6.0")
)

lazy val root = Project(id = "root", base = file("."),
    settings = Seq(
        name := "spark",
        version := "0.1",
        scalaVersion := "2.10.6",

        // Change default location of resources
        //resourceDirectory in Compile := file("./src/main/resources"),
        // Change default location of unmanaged libraries
        //unmanagedBase in Compile := file("./lib"),

        // Apparently sbt mistakes different versions of scala and this fixes it
        fork := true,

        ivyScala := ivyScala.value map {
            _.copy(overrideScalaVersion = true)
        },
        libraryDependencies  ++= Seq(
            "org.scalanlp" %% "breeze" % "0.12",
            "org.scalanlp" %% "breeze-natives" % "0.12",
            "org.scalanlp" %% "breeze-viz" % "0.12",

            "org.apache.spark" %% "spark-core" % versions('spark) withSources() withJavadoc()
        ),

        resolvers ++= Seq(
            "JBoss Repository" at "http://repository.jboss.org/nexus/content/repositories/releases/",
            "Spray Repository" at "http://repo.spray.cc/",
            "Cloudera Repository" at "https://repository.cloudera.com/artifactory/cloudera-repos/",
            "Akka Snapshot Repository" at "http://repo.akka.io/snapshots/",
            "Apache HBase" at "https://repository.apache.org/content/repositories/releases",
            "Twitter Maven Repo" at "http://maven.twttr.com/",
            "scala-tools" at "https://oss.sonatype.org/content/groups/scala-tools",
            "Typesafe repository" at "http://repo.typesafe.com/typesafe/releases/",
            "Second Typesafe repo" at "http://repo.typesafe.com/typesafe/maven-releases/",
            "Mesosphere Public Repository" at "http://downloads.mesosphere.io/maven",
            "snapshots" at "http://oss.sonatype.org/content/repositories/snapshots",
            "Sonatype Releases" at "https://oss.sonatype.org/content/repositories/releases/",
            Resolver.sonatypeRepo("public")
        )
    )
            ++ packSettings
            ++ Seq(
        packMain := Map(
            "template" -> "spark"
        )
    )
)