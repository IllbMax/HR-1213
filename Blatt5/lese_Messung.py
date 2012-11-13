import glob
import re

# Python --version >=  2.7 


def readData(file):
  """
  Liesst alle Messwerttupel (Berechnungszeit, Speed) aus file ein
  und gibt alle die Tupel als Liste aus.
  """
  f = open(file, "r") 
  data = []
  time = []
  speed = []
  for line in f:
    # Bechrenungszeit
    if re.match("Berechnungszeit:\s+[0-9]+.[0-9]+ s" ,line):
      time.append(re.search("(?<= )[0-9]+.[0-9]+", line).group(0))
    # Geschwindigkeit
    if re.match("Speed:\s+[0-9]+.[0-9]+ MFlop/s" ,line):
      speed.append(re.search("(?<= )[0-9]+.[0-9]+", line).group(0))
     
  f.close()
  for i in range(len(time)):
    data.append((time[i], speed[i]))
  return data
  

def parseData(inputfiles = "output/*",
     outputfile = "messtabelle.txt", outputfile_avg = "messtabelle_avg.txt",
     trenneDaten = True):
  """
  Laed alle Daten aus inputfiles, schreibt die Parameter und Daten von readData(file)
  in outputfile, dabei koennen die einzelnen Parameterstaetze durch newline 
  getrennt werden (wenn trenneDaten = True)   
  """
  files = glob.glob(inputfiles)
  files.sort()
  
  tab = open(outputfile, "w")
  tab.write("Threads\tMethod\tInterlines\tinf_func\tTermination\tTerm-Param\tZeit\tSpeed" + "\n")
  
  param = {}
  
  lparam = ""
  for i in files:
    f = re.split("/", i)[-1]
    f = re.split("#", f)[0]
    s = ""
    p = ""
    # Parameter bestimmen
    for d in re.split("_", f):
      p += d + "\t"
    if trenneDaten and lparam != "" and lparam != p:
      tab.write("\n")
    lparam = p
    
    # Zeitbestimmen
    for (time,speed) in readData(i):
      if p not in param:
        param[p] = [(time, speed)]
      else:
        param[p].append((time, speed))
        
      s += p + time + "\t" + speed + "\n"
      
    s = s.strip()
    
    
      
    tab.write(s + "\n")
    
  
  tab.close()
  tab2 = open(outputfile_avg, "w")
  out = []
  for (p,l) in param.iteritems():
    n = len(l)
    x = 0
    xx = 0
    y = 0
    yy = 0
    for (t, s) in l:
      x += float(t)
      xx += float(t)**2
      y += float(s)
      yy += float(s)**2
    
    xmid = x / n
    xvar = ((xx / n) - (x/n)**2) / (n+1) 
    ymid = y / n
    yvar = ((yy / n) - (y/n)**2) / (n+1) 
    
    out.append(p + str(xmid) + "\t" + str(xvar) +
     "\t" + str(ymid) + "\t" + str(yvar))
  
  out.sort()
  tab2.write("Threads\tMethod\tInterlines\tinf_func\tTermination\tTerm-Param\tZeit\tVar(t)\tSpeed\tVar(s)" + "\n")
  for s in out:
    tab2.write(s + "\n")
  
  tab2.close()  
  

#########################
#####     main      #####
#########################    
    
if __name__ == "__main__":
  ###   Steuerparameter   ###
  trenneDaten = True  # bei true zwischen jedem Parametersatzt eine Zeile frei gelassen 
  inputfiles = "output/*"
  outputfile = "messtabelle.txt"
  outputfile_avg = "messtabelle_avg.txt"
  
  parseData(inputfiles, outputfile, outputfile_avg, trenneDaten)
  
  