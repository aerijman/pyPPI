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
-- Table structure for table `perAtomASA`
--

DROP TABLE IF EXISTS `perAtomASA`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `perAtomASA` (
  `perAtomASA_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `Chain` varchar(5) NOT NULL DEFAULT '',
  `Residue` varchar(4) NOT NULL DEFAULT '',
  `ResId` int(11) NOT NULL,
  `Symbol` varchar(4) NOT NULL DEFAULT '',
  `Atom` varchar(2) NOT NULL DEFAULT '',
  `ASA` float NOT NULL,
  `Bfactor` float NOT NULL,
  `Seperated` char(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`perAtomASA_id`),
  UNIQUE KEY `PDB_3` (`PDB`,`Chain`,`Residue`,`ResId`,`Symbol`,`Seperated`),
  KEY `ResId` (`ResId`),
  KEY `atASA` (`ASA`),
  KEY `sepIndx` (`Seperated`)
) ENGINE=InnoDB AUTO_INCREMENT=5008810 DEFAULT CHARSET=latin1 COMMENT='latin1_swedish_ci';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Ndrieding`
--

DROP TABLE IF EXISTS `Ndrieding`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Ndrieding` (
  `fullDrieding_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(4) NOT NULL DEFAULT '',
  `DonorChain` varchar(2) NOT NULL DEFAULT '',
  `DonorResId` int(11) NOT NULL,
  `DonorSymbol` varchar(4) NOT NULL DEFAULT '',
  `AccChain` varchar(2) NOT NULL DEFAULT '',
  `AccResId` int(11) NOT NULL,
  `AccSymbol` varchar(4) NOT NULL DEFAULT '',
  `Energy` float NOT NULL,
  PRIMARY KEY (`fullDrieding_id`),
  UNIQUE KEY `PDB` (`PDB`,`DonorChain`,`DonorResId`,`DonorSymbol`,`AccChain`,`AccResId`,`AccSymbol`,`Energy`)
) ENGINE=InnoDB AUTO_INCREMENT=54429 DEFAULT CHARSET=latin1 COMMENT='latin1_swedish_ci';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `interfaceDist`
--

DROP TABLE IF EXISTS `interfaceDist`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `interfaceDist` (
  `interfaceDist_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `Chains` varchar(6) NOT NULL DEFAULT '',
  `Chain` varchar(4) NOT NULL DEFAULT '',
  `ResId` int(11) NOT NULL,
  `Symbol` varchar(4) NOT NULL DEFAULT '',
  `Atom` varchar(2) NOT NULL DEFAULT '',
  `MinDist` float NOT NULL,
  PRIMARY KEY (`interfaceDist_id`),
  UNIQUE KEY `PDB` (`PDB`,`Chains`,`Chain`,`ResId`,`Symbol`)
) ENGINE=InnoDB AUTO_INCREMENT=51420 DEFAULT CHARSET=latin1 COMMENT='latin1_swedish_ci';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `NinterfaceAtoms`
--

DROP TABLE IF EXISTS `NinterfaceAtoms`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `NinterfaceAtoms` (
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `Chain` varchar(5) NOT NULL DEFAULT '',
  `Residue` varchar(4) NOT NULL DEFAULT '',
  `ResId` int(11) NOT NULL,
  `Symbol` varchar(4) NOT NULL DEFAULT '',
  `atom` varchar(2) NOT NULL DEFAULT '',
  `diffASA` double NOT NULL DEFAULT '0',
  `pk` int(11) NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`pk`),
  UNIQUE KEY `lookup` (`PDB`,`Chain`,`Residue`,`ResId`,`Symbol`)
) ENGINE=InnoDB AUTO_INCREMENT=66131 DEFAULT CHARSET=latin1 COMMENT='latin1_swedish_ci';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `electrostat`
--

DROP TABLE IF EXISTS `electrostat`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `electrostat` (
  `electrostat2_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(4) NOT NULL DEFAULT '',
  `electro` float NOT NULL,
  PRIMARY KEY (`electrostat2_id`),
  UNIQUE KEY `PDB` (`PDB`)
) ENGINE=InnoDB AUTO_INCREMENT=145 DEFAULT CHARSET=latin1 COMMENT='latin1_swedish_ci';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `nRMSD`
--

DROP TABLE IF EXISTS `nRMSD`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `nRMSD` (
  `nRMSD_id` int(11) NOT NULL AUTO_INCREMENT,
  `Complex` varchar(5) NOT NULL DEFAULT '',
  `Unbound` varchar(5) NOT NULL DEFAULT '',
  `Chain` varchar(5) NOT NULL DEFAULT '',
  `UnboundChain` varchar(5) NOT NULL DEFAULT '',
  `RMSD` float NOT NULL,
  `iRMSD` float NOT NULL,
  `Atoms` int(11) NOT NULL,
  `iAtoms` int(11) NOT NULL,
  PRIMARY KEY (`nRMSD_id`),
  UNIQUE KEY `Complex` (`Complex`,`Unbound`,`Chain`,`UnboundChain`)
) ENGINE=InnoDB AUTO_INCREMENT=290 DEFAULT CHARSET=latin1 COMMENT='latin1_swedish_ci';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `interfaceVDW`
--

DROP TABLE IF EXISTS `interfaceVDW`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `interfaceVDW` (
  `interfaceVDW_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `VDV` float NOT NULL,
  `VDVx6` float NOT NULL,
  `ClashV` float NOT NULL,
  `ClashS` float NOT NULL,
  PRIMARY KEY (`interfaceVDW_id`),
  UNIQUE KEY `PDB` (`PDB`)
) ENGINE=InnoDB AUTO_INCREMENT=146 DEFAULT CHARSET=latin1 COMMENT='latin1_swedish_ci';
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

--
-- Table structure for table `interfacePeriphrial`
--

DROP TABLE IF EXISTS `interfacePeriphrial`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `interfacePeriphrial` (
  `interfacePeriphrial_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `Chain` varchar(5) NOT NULL DEFAULT '',
  `ResId` int(11) NOT NULL,
  `Symbol` varchar(5) NOT NULL DEFAULT '',
  `Peri` float NOT NULL,
  `PropPeri` float NOT NULL,
  PRIMARY KEY (`interfacePeriphrial_id`),
  UNIQUE KEY `PDB` (`PDB`,`Chain`,`ResId`,`Symbol`)
) ENGINE=InnoDB AUTO_INCREMENT=388161 DEFAULT CHARSET=latin1 COMMENT='latin1_swedish_ci';
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

--
-- Table structure for table `proteinComplex`
--

DROP TABLE IF EXISTS `proteinComplex`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `proteinComplex` (
  `proteinComplex_id` int(11) NOT NULL AUTO_INCREMENT,
  `PDB` varchar(5) NOT NULL DEFAULT '',
  `UnboundChainA` varchar(6) NOT NULL DEFAULT '',
  `NameA` varchar(55) NOT NULL DEFAULT '',
  `UnboundChainB` varchar(6) NOT NULL DEFAULT '',
  `NameB` varchar(55) NOT NULL DEFAULT '',
  PRIMARY KEY (`proteinComplex_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='latin1_swedish_ci';
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2013-10-05 22:14:00
