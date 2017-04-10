-- MySQL dump 10.13  Distrib 5.5.32, for debian-linux-gnu (i686)
--
-- Host: localhost    Database: ppi
-- ------------------------------------------------------
-- Server version	5.5.32-0ubuntu0.12.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `donors2`
--

DROP TABLE IF EXISTS `donors2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `donors2` (
  `donors2_id` int(11) NOT NULL AUTO_INCREMENT,
  `Residue` varchar(4) NOT NULL DEFAULT '',
  `Symbol` varchar(5) NOT NULL DEFAULT '',
  `Hydrogen` varchar(5) NOT NULL,
  PRIMARY KEY (`donors2_id`),
  UNIQUE KEY `Residue` (`Residue`,`Symbol`,`Hydrogen`)
) ENGINE=InnoDB AUTO_INCREMENT=19 DEFAULT CHARSET=latin1 COMMENT='latin1_swedish_ci';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `donors2`
--

LOCK TABLES `donors2` WRITE;
/*!40000 ALTER TABLE `donors2` DISABLE KEYS */;
INSERT INTO `donors2` VALUES (18,'-','N','H'),(13,'ARG','NE','HE'),(9,'ARG','NH1','1HH1'),(10,'ARG','NH1','2HH1'),(11,'ARG','NH2','1HH2'),(12,'ARG','NH2','2HH2'),(14,'ASN','ND2','1HD2'),(15,'ASN','ND2','2HD2'),(7,'CYS','SG','HG'),(17,'GLN','NE2','1HE2	'),(8,'HIS','NE2','HE2'),(4,'LYS','NZ','1HZ'),(5,'LYS','NZ','2HZ'),(6,'LYS','NZ','3HZ'),(2,'SER','OG','HG'),(3,'THR','OG1','HG1'),(16,'TRP','NE1','HE1'),(1,'TYR','OH','HH');
/*!40000 ALTER TABLE `donors2` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2013-10-18 11:54:32
