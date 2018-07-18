import socket
# import subprocess
import sys

def formatear_mensaje(s):
	host = socket.gethostname()
	return "{0}: {1}".format(host, s)

def mostrar_error(s):
	sys.stderr.write(formatear_mensaje(s) + "\n")

def mostrar_mensaje(s):
	print(formatear_mensaje(s))
	
if __name__ == "__main__":
	if len(sys.argv) == 2:
		mostrar_mensaje(sys.argv[1])
	mostrar_error("ERROR: no se puede dividir por 0!")
	# p = subprocess.Popen("Rscript --version", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	# for line in p.stdout.readlines():
	#	mostrar_mensaje(line)
