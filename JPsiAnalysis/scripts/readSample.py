#!/usr/bin/env python
from ROOT import *
import sys,os
import xml.etree.ElementTree as ET

class Sample :
  def __init__(self, subchild) :
    self.title = subchild.text
    self.Xsec = eval(subchild.attrib['xsec'] )
  def PrintAll(self) :
    print self.title, self.Xsec

class HistInfo :
  def __init__(self, child) :
    self.SampleList=[]
    self.color = child.attrib['color']
    self.title = child.attrib['title']
    for subchild in child :
      self.SampleList.append(Sample( subchild) )
  def PrintAll(self) :
    print self.color, self.title
    for sample in self.SampleList :
      sample.PrintAll()   
 
class readSample :
  def __init__(self, inputXML) :
    self.HistInfoList =[]
    eltree = ET.parse(inputXML)
    root = eltree.getroot()
    for child in root :
      self.HistInfoList.append( HistInfo ( child ))
  def PrintAll(self) :
    for hinfo in self.HistInfoList :
      hinfo.PrintAll()
  
    
if __name__ == '__main__' : 
  if ( len(sys.argv) != 2 ) :
    os.system(-1)
  rs = readSample(sys.argv[1])
  rs.PrintAll() 
   
