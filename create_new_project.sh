#!/bin/bash

## creat new project
wkdir=$1

## creaste working directory
mkdir -p $wkdir

## data: final report, sample table, centromere position, and duplicate pairs [optional]
## put in this directory
mkdir -p ${wkdir}/data

## 01_initial_call
cp -ru ./01_initial_call $wkdir
mkdir -p ${wkdir}/01_initial_call/finalreport_to_matrix_LRR_and_BAF/RDS

mkdir -p ${wkdir}/01_initial_call/run_iPattern/data
mkdir -p ${wkdir}/01_initial_call/run_iPattern/data_aux
mkdir -p ${wkdir}/01_initial_call/run_iPattern/results

mkdir -p ${wkdir}/01_initial_call/run_PennCNV/data
mkdir -p ${wkdir}/01_initial_call/run_PennCNV/data_aux
mkdir -p ${wkdir}/01_initial_call/run_PennCNV/results

mkdir -p ${wkdir}/01_initial_call/run_QuantiSNP/data
mkdir -p ${wkdir}/01_initial_call/run_QuantiSNP/data_aux
mkdir -p ${wkdir}/01_initial_call/run_QuantiSNP/results

## 02_batch_effect
cp -ru ./02_batch_effect $wkdir

## 03_create_CNVR      
cp -ru ./03_create_CNVR $wkdir

## 04_CNV_genotype
cp -ru ./04_CNV_genotype $wkdir
mkdir -p ${wkdir}/04_CNV_genotype/data
mkdir -p ${wkdir}/04_CNV_genotype/results

## 05_boundary_refinement 
cp -ru ./05_boundary_refinement $wkdir
mkdir -p ${wkdir}/05_boundary_refinement/data
mkdir -p ${wkdir}/05_boundary_refinement/results

## 05a_regenotype_after_refinement
mkdir -p ${wkdir}/05a_regenotype_after_refinement
mkdir -p ${wkdir}/05a_regenotype_after_refinement/data
mkdir -p ${wkdir}/05a_regenotype_after_refinement/results

## 06_performance_assessment
cp -ru ./06_performance_assessment $wkdir

echo "New project directory has been created at: $wkdir"
echo "Please put (or create symbolic link to) input data in the directory: $wkdir/data"

