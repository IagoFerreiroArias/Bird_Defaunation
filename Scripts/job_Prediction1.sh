#!/bin/sh
#SBATCH -N 1 #(solicita un nodo)
#SBATCH -n 1 --mem-per-cpu=64G #con 8 G da fallo
#SBATCH -t 24:00:00 #10 min
#SBATCH -c 30 # cores per task
#SBATCH --mail-type=begin #Envía un correo cuando el trabajo inicia
#SBATCH --mail-type=end #Envía un correo cuando el trabajo finaliza
#SBATCH --mail-type=fail #Envía un correo cuando el trabajo falla
#SBATCH --mail-user=iago.ferreiro.arias@gmail.com #Dirección a la que se envía
R --no-save < Predictions1.R