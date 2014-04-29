# =========================================================================
#  Program:   MATITK: Extending MATLAB with ITK
#  Version:   2.4.02
#  Language:  C++
#
#  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
# Lab, Simon Fraser University. All rights reserved. 
#  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
#  details.
#
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.  See the above copyright notices for more information.
# =========================================================================

$MS='///////////////begin core filter code///////////////////////';
$ME='///////////////end core filter code///////////////////////';

#Step 1 - Read all files, extract interesting regions
#       - Find out what functions are there
opendir(DIR,".");
open(OH,">step1.txt");

while ($f=readdir(DIR)) {
        if ($f=~/(.+)\.cxx$/){
          open(FH,"<$f");
          $intersted=0;   #we are interested in code surrounded by BeginCodeSnippet and EndCodeSnippet
	  print OH "\nvoid $1(){\n";
          push @fcns,$1;
	  print OH "$MS\n";
	  $filtername='filter';
          while (<FH>){
		if(/BeginCodeSnippet/) {$interested=1;}
		elsif(/EndCodeSnippet/) {$interested=0;}
		elsif ($interested) {
                  print OH;
                  /^(\s+)/; 
                  $indent=$1;}
		if ($interested && /writer\-\>SetInput\s*\(\s*([A-Za-z0-9]*)/){  #find out what the main filter's name is
	  		$filtername=$1; 
        	}
	        
          }
          print OH $indent,"pixelContainer = $filtername->GetOutput()->GetPixelContainer();\n$ME\n}\n";
          close FH;
        } #endif cxx     
}#end while

my $disp_temp='const char* OPCODE[]={"REPLACEME1"};

const char* OPNAME[]={"REPLACEME2"};

pt2Function  OPFCNARRAY[]={&REPLACEME3};';
my $insert1=join '","',@fcns;
$insert1=~s/[a-z]//g;
my $insert2=join "\",\n\"",@fcns;
my $insert3=join ', &',@fcns;
$disp_temp=~s/REPLACEME1/$insert1/;
$disp_temp=~s/REPLACEME2/$insert2/;
$disp_temp=~s/REPLACEME3/$insert3/;
print OH "\/*\n\/\/Insert to dispatcher...\n$disp_temp\n";

closedir(DIR);
close OH;

#Step 2 - Remove redundant and rearrange includes
open(FH,"<step1.txt");
while (<FH>){
    my $line=$_;
    $step2flag=0;
    @hated=('InputImageType;','OutputImageType;','InputPixelType','OutputPixelType','writer->SetInput','writer->Update','ReaderType','#include');
    foreach $array_element(@hated)
    {
       if ($line=~/$array_element/) {$step2flag=1;}  #ignore the line if it's one of the hated strings
    }
    unless ($step2flag) {$accum.=$line;}
    if ($line=~/#include/) {   #save all the includes and put to the top of file
 	unless ($includes=~/$line/) {
		$includes=$includes.$line;}
	}
     }
open(OH,">step2.txt");
print OH $includes;
print OH $accum;
close OH;
close FH;

#step 3 - Replace reader->GetOutput() with importFilter[IMPORTFILTERA]->GetOutput() 
#         Replace argv with constant names

open(FH,"<step2.txt");
open(OH,">step3.txt");
while (<FH>){
        if (/$MS/) {$count=0;}
	s/reader\-\>GetOutput\(\)/importFilter\[IMPORTFILTERA\]\-\>GetOutput\(\)/;
        if (/argv/ && /->/){   #replace argv with constant
           /([^-]*->)([^\(]*)/;
	   print OH '//',$_;
	   my $consname=uc($2);
           print OH "$1$2\($consname\);\n";
           $count++;
	}
	elsif (/argv/){
	   print OH '//',$_;
	   s/atoi/\(unsigned int\)/;
           s/atof//;
           /^\s+([a-zA-Z0-9]+)/;
           my $consname=uc($1);
	   s/argv\[\d+\]/PARAM$count$consname/;
           print OH "$_";
           $count++;
	}
        else {
	   print OH;
	}
}
close OH;
close FH;


#step 4 - Create parameters container
open(FH,"<step3.txt");
open(OH,">step4.txt");

my $template='	const char* PARAM[]={"REPLACEME1"};
	const char* SUGGESTVALUE[]={REPLACEME2};
	const int nParam = sizeof(PARAM)/sizeof(*PARAM);
	ParameterContainer paramIterator(PARAM,SUGGESTVALUE,nParam);';

while (<FH>){
   if (/$MS/) {
      @constants=(); $accum='';$ROI=1;}
   elsif (/$ME/){

      my $insert=join '","',@constants;
      $_=$template;
      s/REPLACEME1/$insert/;
      my $sugvalue='"",'x@constants;
      $sugvalue=~s/,$//;
      s/REPLACEME2/$sugvalue/;
      print OH "$_\n";
      for ($count=0; $count<@constants; $count++)
	{
	  my $curr=$constants[$count];
	  print OH "double $curr=paramIterator.getCurrentParam($count);\n";
      }
      print OH "$accum";
      $ROI=0;
      print OH "$ME\n";
      next;
   }
   elsif (/\s+filter(.*?)Set(.*?)\(\s*([A-Za-z0-9]+)\s*\)/){
      push @constants,$3;
   }
   elsif (/\(\s*([A-Z]+[0-9]*[A-Z]*)\s*\)/){
      push @constants,$1;
   }
   if ($ROI) {
      $accum.=$_;
   }
   else {print OH;}
}
close OH;
close FH;