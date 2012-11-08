import glob
import re


def readData(file):
  f = open(file, "r") 
  data = []
  for line in f:
    # Bechrenungszeit
    if re.match("Berechnungszeit:\s+[0-9]+.[0-9]+ s" ,line):
      data.append(re.search("(?<= )[0-9]+.[0-9]+", line).group(0))
    # Geschwindigkeit
    if re.match("Speed:\s+[0-9]+.[0-9]+ MFlop/s" ,line):
      data.append(re.search("(?<= )[0-9]+.[0-9]+", line).group(0))
      
  f.close()
  return data
  



dir = "output/*"

files = glob.glob(dir)


tab = open("messtabelle.txt", "w")
tab.write("Threads\tMethod\tInterlines\tinf_func\tTermination\tTerm-Param\tZeit\tSpeed" + "\n")

for i in files:
  f = re.split("/", i)[-1]
  f = re.split("#", f)[0]
  s = ""
  # Parameter bestimmen
  for d in re.split("_", f):
    s += d + "\t"
  
  # Zeitbestimmen
  for d in readData(i):
    s += d + "\t"
  s = s.strip()
    
  tab.write(s + "\n")
  

tab.close()