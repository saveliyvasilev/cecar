#!/usr/bin/Rscript

computeData <- function(){
	cat(sprintf("Hola, soy un script en R. Voy a intentar cargar un input de %s.\n", input_file))
	data <- read.csv(input_file)
	print(data)
	print("Maximo salario")
	salmax <- max(data$salary)
	print(salmax)
	print("Gente de IT")
	retval <- subset(data, dept=="IT")
	print(retval)
	cat(sprintf("Trato de guardar en %s.\n", output_file))
	write.csv(retval, output_file)
}

show_random <- function(){
	x1 <- runif(1, 0, 100)
	sprintf("Numero aleatorio: %f", x1)
	print(input_file)
}


# Leo los argumentos de programa, es necesario para luego
# agregar los resultados. Se pasa el input file y outputfile.

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2){
	stop("Faltan los archivos de input y output. Debe ser: $ Rscript <ejecutable> <input> <output>")
} else {
	input_file <- args[1]
	output_file <- args[2]
	cat(sprintf("INPUT: %s\n", input_file))
	cat(sprintf("OUTPUT: %s\n", output_file))
}

computeData()
show_random()


