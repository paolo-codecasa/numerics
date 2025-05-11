# Datei, die kompiliert werden soll (ohne Endung)
BINARY := epidemy

# Header Datei, in der alle notwendigen Numerik Funktionen sind (../ = einen Ordner hoch gehen)
HEADER := ../../my_numerics

# Argumente, mit dem das Programm ausgeführt wird,
# falls das Programm keine Argumente aus der Konsole übernimmt, kann das leer gelassen werden
ARGS   :=

# Der zu verwendende C Compiler
CC     := gcc

# Flags, die für das Kompilieren benötigt werden
CFLAGS := -O3 -lm -lgsl -lgslcblas -l lapacke

# Flags, die zum Debugging benötigt werden
DFLAGS := -Wall -Wextra -Wshadow -pedantic  -fstack-protector

# Flags, die für das linken von .o Dateien benötigt werden
LFLAGS := -flto

# -lm                -> math.h / tgmath.h
# -lgsl / -lgslcbla  -> gsl.h
# -l lapacke         -> lapacke.h
#-----------------------
# -O3 				 -> Optimiert den Code, es gibt -O0 < -O1 < -O2 < -O3 < -Ofast,
# 						wobei bei mehr Optimierungen (in Richtung -Ofast) mehr Präzision verloren geht
# -Wall				 -> Aktiviert Compiler-Warnflags
# -flto				 -> Optimierung der Verlinkung von .o Dateien

# Struktur:
#-----------------------
# ZIEL: ABHÄNGIGKEIT
#	BEFEHLE
#-----------------------
# Benutzung:
#-----------------------
# make  	-> Nur kompilieren
# make run 	-> Kompilieren und ausführen
# make plot -> Ausführen von Pythonskripten
# make clean-> Nicht notwendige Dateien löschen

# Nachfolgend werden die oben definierten Variablen benutzt, um das Kompilieren zu vereinfachen

# Kompilierung der Zieldatei, Verlinkung mit der Headerdatei
$(BINARY): $(BINARY).c $(HEADER).c $(HEADER).o
	gcc $(LFLAGS) -o $(BINARY) $(BINARY).c $(HEADER).o $(DFLAGS) $(CFLAGS)

# Kompilierung der Headerdatei
$(HEADER).o: $(HEADER).c
	gcc -c $(HEADER).c -o $(HEADER).o $(DFLAGS) $(CFLAGS)

# Wenn das Programm existiert, führe es aus, wenn nicht kompiliere und führe es dann aus
run: $(BINARY)
	./$< $(ARGS)

# Führe Programm aus und danach sofort das angegebene Python Script
plot: run
	python3 plot.py

# Entfernen von temporären und nicht notwendigen Dateien
clean:
	rm -f *.o
	rm -f $(HEADER).o
	rm -f *.exe
	rm -f $(BINARY)
	rm -f *.dat *.csv *.txt