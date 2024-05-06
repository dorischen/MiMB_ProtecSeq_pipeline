/**
 * @(#)SAMManipulator.java
 *
 * input: SAM alignment file, read ids
 * output: SAM alignment file with or without specified reads

 *
 * @doris chen
 * @version 170206: addition of RG per sample (ccs)
 * @version 170801: addition of RG per sample (subreads)
 */

import java.io.*;
import java.util.*;                    // for ArrayList
import java.util.regex.*;              // for regular expression methods (s. also RegexTestHarness.java)


public class SAMManipulator
{  static File inputFileS, inputFileI;
   static String option, outfileBase, format, library, machine;  
   static Integer readCount, rgCount;
  
   public static void main (String args[]) throws Exception
   {
      // default
      option="";
      outfileBase="";
      format="";
      library="";
      readCount=0;
      machine="";
      rgCount=0;
    
      try
      {  if (args.length>=4)
			      {  Pattern pSam = Pattern.compile("^-sam");
			         Pattern pIds = Pattern.compile("^-ids");
			         Pattern pOption = Pattern.compile("^-option");
			         Pattern pOut = Pattern.compile("^-out");
			         Pattern pFormat = Pattern.compile("^-format");
			         Pattern pLib = Pattern.compile("^-lib");
			         Pattern pCount = Pattern.compile("^-count");
			         Pattern pMach = Pattern.compile("^-mach");
			         Pattern pRGCount = Pattern.compile("^-rgcount");
			       
            for (int i=0; i<args.length; i++)
			         {  Matcher mSam = pSam.matcher(args[i]);
			            if (mSam.find())
			            {  inputFileS = new File(args[i+1]);
			               if ( !inputFileS.exists() ) throw new Exception();
			            }
               else
               {  Matcher mIds = pIds.matcher(args[i]);
						            if (mIds.find())
						            {   inputFileI = new File(args[i+1]);
                      if ( !inputFileI.exists() ) throw new Exception();
						            }
						            else
						            {  Matcher mOption = pOption.matcher(args[i]);
									            if (mOption.find())
									            {  option = args[i+1];
									               if ( !option.equals("delete") && !option.equals("keep") && !option.equals("addRG") && !option.equals("addRGH") && !option.equals("addRGHS") && !option.equals("subsample")) throw new Exception();
									            }  else
									            {  Matcher mOut = pOut.matcher(args[i]);
												            if (mOut.find())
												            {  outfileBase = args[i+1];															               
												            } else
												            {  Matcher mFormat = pFormat.matcher(args[i]);
															            if (mFormat.find())
															            {  format = args[i+1];
									                     if ( !format.equals("ccs") && !format.equals("subread") && !format.equals("read") && !format.equals("pe") && !format.equals("se") ) throw new Exception();															               
															            }  else
															            {  Matcher mLib = pLib.matcher(args[i]);
																		            if (mLib.find())
																		            {  library = args[i+1];												                    												               
																		            }  else
																		            {  Matcher mCount = pCount.matcher(args[i]);
																					            if (mCount.find())
																					            {  readCount = Integer.parseInt(args[i+1]);												                    												               
																					            } else																				            
																					            {  Matcher mMach = pMach.matcher(args[i]);
																								            if (mMach.find())
																								            {  machine = args[i+1];												                    												               
																								            } else
																								            {  Matcher mRGCount = pRGCount.matcher(args[i]);
																											            if (mRGCount.find())
																											            {  rgCount = Integer.parseInt(args[i+1]);												                    												               
																											            } 
																								            }
																					            }
																					         }
																		         }
															         }
												         }
									         }
															}
			         }
            
            if(option.equals("addRG") || option.equals("addRGH") || option.equals("addRGHS"))
            	 add_rg_to_sam(inputFileS, outfileBase, format, library, machine, rgCount);
            else if(option.equals("subsample"))
            	 subsample(inputFileS, readCount, outfileBase);
            else
              process_sam(inputFileS, inputFileI, option, outfileBase, format);

			      }
			      else throw new Exception();
			   }
			   catch (Exception exc)
			   {  System.out.println("");
			      System.out.println(">java SAMManipulator -sam SAM_FILE -option delete|keep|addRG|addRGH|addRGHS|subsample -out OUTFILE_BASE [-ids ID_FILE] [-format ccs|subread|read|pe|se -lib LIBRARY_NAME -mach SEQUENCING_MACHINE(e.g. HISEQ2500 or SEQUEL) -rgcount NR_OF_READGROUPS] [-count READ(PAIR)_COUNT]");
			      System.out.println("ad ID_FILE: read (or hole) ids, one line each");
			      System.out.println("ad option .. delete -> reads deleted from sam; keep -> reads kept (rest of reads deleted); addRG .. add RG in header and in alignments; addRGH .. add RG in header only; addRGHS .. add RG in header only, 'flowcell' information from first substring in read name separated by dot");
			      System.out.println("ad OUTFILE_BASE: file name for resulting sam file (.sam appended)");
			      System.out.println("ad format: 'ccs' .. in case of ccs alignment, extraction by hole ids (optional); 'subread' not implemented yet; read .. ids taken as such; pe .. PAIRED in RG header");
			      System.out.println("ad LIBRARY_NAME: added to RG header line, e.g. 'SK1-BY_msh2d_ccs2'");
			      System.out.println("For extraction of alignments per read:");
			      System.out.println("java SAMManipulator -sam SAM_FILE -option delete|keep -format ccs|subread -out OUTFILE_BASE -ids ID_FILE");			  
			      System.out.println("For addition of read group information:");
			      System.out.println("java SAMManipulator -sam SAM_FILE -option addRG|addRGH|addRGHS -format ccs|subread|read|pe|se -lib LIBRARY_NAME -mach SEQUENCING_MACHINE -rgcount NR_OF_READGROUPS -out OUTFILE_BASE");			   
			      System.out.println("For subsampling certain number of reads or read pairs randomly:");
			      System.out.println("java SAMManipulator -sam SAM_FILE -option subsample -count READ(PAIR)_COUNT -out OUTFILE_BASE");			     
			     
         System.out.println("Error: " + exc.toString() + '\n');
			      System.out.println("");
			   }

   }
   
   public static void subsample(File inputFileS, Integer readCount, String outfileBase)
	  {  try
			   {  ReadFile rfS = new ReadFile(inputFileS);
         Pattern pTab = Pattern.compile("\t"); 
         Pattern pSlash = Pattern.compile("/"); 
         Pattern pAt = Pattern.compile("^@");
       
			      // for output file 
         File resultFile = new File(rfS.getCanonParent() + "/" + outfileBase +"_sub"+readCount+".sam");
         PrintWriter outS = new PrintWriter(new FileWriter(resultFile));
         
         // go through sam file
         BufferedReader inS = new BufferedReader(new FileReader(inputFileS));
         String currentLineS = inS.readLine();
         // header lines
         while(currentLineS!=null)
         {  Matcher mAt = pAt.matcher(currentLineS);
            if(mAt.find())
            {  outS.println(currentLineS);
            	  currentLineS = inS.readLine();
            }
            else
            	 break;
         }
         
         ArrayList<String> rowList = new ArrayList<String>(1000000);
         while(currentLineS!=null)
         {  rowList.add(currentLineS);
            currentLineS = inS.readLine();
			      }
         inS.close();
         
         Collections.shuffle(rowList);
         
         // print random rows
         int rowCounter=0;
         for(String row: rowList)
         {  if(rowCounter<=readCount)
            {  outS.println(row);   
            	  rowCounter++;
            } else
            	 break;      	
         }
         outS.close();
         
         System.out.println(rowCounter+" random alignment(s) selected.");
         System.out.println("");

     }
			  catch (Exception exc)
			  { System.out.print("Error: " + exc);
			    System.exit(1);
			  }
		}
		
 
   public static void process_sam(File inputFileS, File inputFileI, String option, String outfileBase, String format)
	  {  try
			   {  boolean delete=false;
			      String suffix = "retained";
			   	  if(option.equals("delete"))
			   	  {  delete = true;
			   	     suffix = "deleted";
			   	  }
			   	  boolean ccs=false;
			   	  if(format.equals("ccs"))
			   	    ccs = true;
			   	  			   	  
			   	  // get reads
         ReadFile rfI = new ReadFile(inputFileI);
         BufferedReader inI = new BufferedReader(new FileReader(inputFileI));
         
         ArrayList<String> idList = new ArrayList(10000);
         String currentLineI = inI.readLine();
         if(ccs)
         {  while(currentLineI!=null)
			         {  idList.add(currentLineI.trim().split("/ccs")[0]);
			         	  currentLineI = inI.readLine();
			         }
         } else
         {  while(currentLineI!=null)
			         {  idList.add(currentLineI.trim());
			         	  currentLineI = inI.readLine();
			         }
         }
         inI.close();
         System.out.println(inputFileI + " read.");

         // go through sam file
         ReadFile rfS = new ReadFile(inputFileS);
         Pattern pTab = Pattern.compile("\t"); 
         Pattern pSlash = Pattern.compile("/"); 
         Pattern pAt = Pattern.compile("^@");
       
			      // for output file 
         File resultFile = new File(rfS.getCanonParent() + "/" + outfileBase +".sam");
         PrintWriter outS = new PrintWriter(new FileWriter(resultFile));
         
         BufferedReader inS = new BufferedReader(new FileReader(inputFileS));
         String currentLineS = inS.readLine();
         // header lines
         while(currentLineS!=null)
         {  Matcher mAt = pAt.matcher(currentLineS);
            if(mAt.find())
            {  outS.println(currentLineS);
            	  currentLineS = inS.readLine();
            }
            else
            	 break;
         }
         
         String currentRead="";
         int rowCounter=0;
         if(delete)
         {  while(currentLineS!=null)
			         {  // go through rows, split into items
			            String[] row = pTab.split(currentLineS);
			            currentRead = row[0];
			            if(!idList.contains(currentRead))
			            {  outS.println(currentLineS.trim()); 
			            	  rowCounter++;
			            }
			            currentLineS = inS.readLine();
						      }
         } else
         {  if(ccs)
            {  while(currentLineS!=null)
						         {  // go through rows, split into items
						            String[] row = pTab.split(currentLineS);
						            String[] name = pSlash.split(row[0]);
						            currentRead = name[1];
						            if(idList.contains(currentRead))
						            {  outS.println(currentLineS.trim()); 		
						            	  rowCounter++;    			            	             	      	
						            } 
						            currentLineS = inS.readLine();
									      }
            } else
            {  while(currentLineS!=null)
						         {  // go through rows, split into items
						            String[] row = pTab.split(currentLineS);
						            currentRead = row[0];						            
						            if(idList.contains(currentRead))
						            {  outS.println(currentLineS.trim()); 		
						            	  rowCounter++;    			            	             	      	
						            } 
						            currentLineS = inS.readLine();
									      }
            }
         }
         inS.close();
         outS.close();
         
         System.out.println(idList.size()+" read(s) "+suffix+".");
         System.out.println(rowCounter+" alignment(s) found.");
         System.out.println("");

     }
			  catch (Exception exc)
			  { System.out.print("Error: " + exc);
			    System.exit(1);
			  }
		}
		
		public static void add_rg_to_sam(File inputFileS, String outfileBase, String format, String library, String machine, Integer rgCount)
	 {  try
			  {  String suffix = "_rg";
			   	 boolean ccs=false, subread=false;
			   	 if(format.equals("ccs"))
			   	 {  ccs = true;
			   	 } else
			   	 if(format.equals("subread"))
			   	 {  subread = true;
			   	 }
			   	 
			   	 // patterns for parsing
        Pattern pTab = Pattern.compile("\t"); 
        Pattern pSlash = Pattern.compile("/"); 
        Pattern pSpace = Pattern.compile(" ");   
        Pattern pAt = Pattern.compile("^@");   
        Pattern pRG = Pattern.compile("^@RG");   
        Pattern pRG2 = Pattern.compile("RG:Z:");   
        Pattern pPG = Pattern.compile("^@PG");
        Pattern pCL = Pattern.compile("CL:");
        Pattern pUnder = Pattern.compile("_");
        Pattern pSA = Pattern.compile("SA:Z:");
        Pattern pHyphen = Pattern.compile(":");
        Pattern pDot = Pattern.compile("\\.");
			   	
			   	 // load sam file
        ReadFile rfS = new ReadFile(inputFileS);
        BufferedReader inS0 = new BufferedReader(new FileReader(inputFileS));
        String currentLineS = inS0.readLine();
        HashSet<String> fcrgList = new HashSet(10);
			     HashSet<String> idList = new HashSet(10000);
        if(option.equals("addRGH"))
        {  // get unique RGs and flowcell ids
           String fc="";
           while(currentLineS!=null)
			        {  Matcher mAt = pAt.matcher(currentLineS);
			           if(!mAt.find())
			           {   // get flowcell id
			           	   String[] lineSplit = pTab.split(currentLineS.trim());
			               String[] fcSplit = pHyphen.split(lineSplit[0]);
			               fc = fcSplit[0]+"_"+fcSplit[1]+"_"+fcSplit[2];
			               // get RG
			               String[] rgSplit = pRG2.split(currentLineS.trim());
			               fcrgList.add(fc+"/"+rgSplit[1].trim());  // RG tag expected at end of line (if not, add split at tabs) 			             
			           }
					         currentLineS = inS0.readLine();			        
			        }
			        inS0.close();
        } else
        if(option.equals("addRGHS"))
        {  // get unique RGs and flowcell ids
           String fc="";         
           while(currentLineS!=null)
			        {  Matcher mAt = pAt.matcher(currentLineS);			          
			           if(!mAt.find() && fcrgList.size()<rgCount)
			           {   // get flowcell id			               
			           	   String[] lineSplit = pTab.split(currentLineS.trim());			           	   
			               String[] fcSplit = pDot.split(lineSplit[0]);			              
			               fc = fcSplit[0];			              
			               // get RG			              
			               String[] rgSplit = pRG2.split(currentLineS.trim());
			               fcrgList.add(fc+"/"+rgSplit[1].trim());  // RG tag expected at end of line (if not, add split at tabs) 				                			                      
			           }
					         currentLineS = inS0.readLine();			        
			        }
			        inS0.close();
        } else
        {  // get unique hole ids of aligned reads first
			        String idSource="";
			        while(currentLineS!=null)
			        {  Matcher mAt = pAt.matcher(currentLineS);
			           if(!mAt.find())
			           {   String[] lineSplit = pTab.split(currentLineS.trim());
			               if(!lineSplit[1].equals("4"))
			               {  if(ccs)
			                  {  idSource = lineSplit[0].replaceAll("/","_");
			                  	  idList.add(idSource);
			                  } else
			                  if(subread)
			                  {  idSource = lineSplit[0].replaceAll("/","_");
			                  	  String[] idSplit = pUnder.split(idSource);
			                  	  idList.add(idSplit[0]+"_"+idSplit[1]+"_"+idSplit[2]+"_"+idSplit[3]+((lineSplit[1].equals("0") | lineSplit[1].equals("2048")) ? "_sr_fwd" : "_sr_rev"));
			                  }
			               }
			           }
					         currentLineS = inS0.readLine();			        
			        }
			        inS0.close();
        }
       
        // for output file 
        String resultFileName="";
        resultFileName = rfS.getCanonParent() + "/" + outfileBase + suffix + ".sam";
        File resultFile = new File(resultFileName);
        PrintWriter outS = new PrintWriter(new FileWriter(resultFile));       
               
        // go through sam header lines
        String aligner="", pl="", pm="", idString;
        boolean plusRG=false, plusAligner=false;
        if(format.equals("ccs") || format.equals("subread"))
        {  pl = "PACBIO";                   	 
        } else
        {  pl = "ILLUMINA";
        }  
        pm = machine;
        
        BufferedReader inS = new BufferedReader(new FileReader(inputFileS));
        currentLineS = inS.readLine();
        while(currentLineS!=null)
        {  Matcher mAt = pAt.matcher(currentLineS);
           if(mAt.find())
           {  Matcher mRG = pRG.matcher(currentLineS);
		            if(!mRG.find())
           	  {  Matcher mPG = pPG.matcher(currentLineS);
					            if(mPG.find())
					            {  if(!plusRG)
					            	  {  if(option.equals("addRGH") || option.equals("addRGHS"))
					            	  	  {  // print new RG lines (one per read group)
																	         for(String fcrg: fcrgList)
																	         {  String[] fcrgSplit = pSlash.split(fcrg);
																	            outS.println("@RG"+"\t"+"ID:"+fcrgSplit[1]+"\t"+"DS:READTYPE="+format.toUpperCase()+"\t"+"LB:"+library+"\t"+"PL:"+pl+"\t"+"PM:"+pm+"\t"+"PU:"+fcrgSplit[0]+"\t"+"SM:"+fcrgSplit[1]);
																	            plusRG = true;
																	         }			
					            	  	  } else
					            	  	  {  // print new RG lines (one per read)
					            	  	     for(String id: idList)
																	         {  String[] idSplit = pUnder.split(id);
																	            outS.println("@RG"+"\t"+"ID:"+(idSplit[3]+"_"+idSplit[4]+"_"+idSplit[5])+"\t"+"DS:READTYPE="+format.toUpperCase()+"\t"+"LB:"+library+"\t"+"PL:"+pl+"\t"+"PM:"+pm+"\t"+"PU:"+(idSplit[0]+"_"+idSplit[1]+"_"+idSplit[2])+"\t"+"SM:"+(idSplit[3]+"_"+idSplit[4]+"_"+idSplit[5]));
																	            plusRG = true;
																	         }						            	  	     
					            	  	  }              
					            	  } 
					            	  if(!plusAligner)
					               {  // get alignment tool name
			           	        String[] pgRow = pCL.split(currentLineS);
								               String[] pgRow2 = pSpace.split(pgRow[1]);
								               aligner = pgRow2[0];
								               System.out.println(aligner);
								               plusAligner = true;
					               }
  			            }	
					            outS.println(currentLineS);  // print all header lines except RG lines like in original file
		            } else  // in case of existing RG lines
		            {  Matcher mPG = pPG.matcher(currentLineS);
					            if(mPG.find())
					            {  if(!plusRG)
					            	  {  if(option.equals("addRGH") || option.equals("addRGHS"))
					            	  	  {  // print new RG lines (one per read group)
					            	  	     int fCounter=0;
																	         for(String fcrg: fcrgList)
																	         {  String[] fcrgSplit = pSlash.split(fcrg);
																	            outS.println("@RG"+"\t"+"ID:"+fcrgSplit[1]+"\t"+"DS:READTYPE="+format.toUpperCase()+"\t"+"LB:"+library+"\t"+"PL:"+pl+"\t"+"PM:"+pm+"\t"+"PU:"+fcrgSplit[0]+"\t"+"SM:"+fcrgSplit[1]);
																	            plusRG = true;
																	            fCounter++;
																	            if(fCounter<rgCount)
																	              currentLineS=inS.readLine();  // jump to next line
																	         }			
					            	  	  } else
					            	  	  {  // print new RG lines (one per read)
					            	  	     for(String id: idList)
																	         {  String[] idSplit = pUnder.split(id);
																	            outS.println("@RG"+"\t"+"ID:"+(idSplit[3]+"_"+idSplit[4]+"_"+idSplit[5])+"\t"+"DS:READTYPE="+format.toUpperCase()+"\t"+"LB:"+library+"\t"+"PL:"+pl+"\t"+"PM:"+pm+"\t"+"PU:"+(idSplit[0]+"_"+idSplit[1]+"_"+idSplit[2])+"\t"+"SM:"+(idSplit[3]+"_"+idSplit[4]+"_"+idSplit[5]));
																	            plusRG = true;
																	         }						            	  	     
					            	  	  }              
					            	  } 
					            }
		            }  	  	
           	  currentLineS = inS.readLine();
           }
           else
           	 break;
        }
       
        // add tag to all alignment rows 
        String currentRead="", rgSplit="", saSplit="", outLine="";
        int rowCounter=0;
        
        // look for already existing RG tag
        int rgIndex=0, saIndex=19;
        boolean rg=false, sa=false;
        Matcher mRG2 = pRG2.matcher(currentLineS);
        if(mRG2.find())
		      {  String[] currentLineSplit = pTab.split(currentLineS);
		         for(int i=0; i<currentLineSplit.length; i++)
		         {   Matcher mRG2b = pRG2.matcher(currentLineSplit[i]);
		             if(mRG2b.find())
		             {	 rgIndex = i;
		                rg = true;			               
		             	  break;	
		             }		             
		         }
		      }
		      
		      if(ccs)
		      {  while(currentLineS!=null)
					      {  outLine = currentLineS.trim();
					         String[] currLineSplit = pTab.split(outLine);
					         
					         Matcher mSA = pSA.matcher(outLine);
		            if(mSA.find())
		            {  saSplit = currLineSplit[saIndex];  // should be added at end of line (otherwise problems with loading into R)
		            	  outLine = outLine.replace((saSplit+"\t"),"");		// print line without SA tag
		            	  sa = true;
		            }
		            
					         if(rg)
					         {  rgSplit = currLineSplit[rgIndex];
					            outLine = outLine.replace((rgSplit+"\t"),"");		// print line without RG tag
					         } 
					         
					         String[] nameSource = pUnder.split(currLineSplit[0].replaceAll("/","_"));
					         currentRead = nameSource[3]+"_"+nameSource[4]+"_"+nameSource[5];
					         outLine = outLine+"\t"+"RG:Z:"+currentRead; 
					         	
					         if(sa)
					         {  outLine = outLine+"\t"+saSplit;
					         		 sa = false;
					         }
					         
					         outS.println(outLine);					         	
					         rowCounter++;    
					         				            	             	      	
		            currentLineS = inS.readLine();						     
		         }        	
		      } else
		      if(subread)
		      {  while(currentLineS!=null)
		      	  {  String[] currLineSplit = pTab.split(currentLineS);
		      	     if(!currLineSplit[1].equals("4"))
		      	     {  if(rg)
								         {  rgSplit = currLineSplit[rgIndex];
								            outS.print(currentLineS.trim().replace((rgSplit+"\t"),""));		// print line without RG tag
								         } else
								         {	 outS.print(currentLineS.trim());  // print whole line and add RG tag
								         }
								         String[] nameSource = pUnder.split(currLineSplit[0].replaceAll("/","_"));
								         currentRead = nameSource[3]+((currLineSplit[1].equals("0") | currLineSplit[1].equals("2048")) ? "_sr_fwd" : "_sr_rev");
								         outS.println("\t"+"RG:Z:"+currentRead); 
								         rowCounter++;   
		      	     } 			            	             	      	
		            currentLineS = inS.readLine();						     
		         }      
		      } else
		      {  while(currentLineS!=null)
		      	  {  outS.println(currentLineS.trim()); 								      	      	     		            	             	      	
		            currentLineS = inS.readLine();						     
		         }   
		         rowCounter = fcrgList.size();   
		      } 
        inS.close();
        outS.close();
        
        System.out.println(rowCounter+" RG tags added.");
        System.out.println("Updated alignment file saved as "+resultFileName);
        System.out.println("");

     }
			  catch (Exception exc)
			  { System.out.print("Error: " + exc);
			    System.exit(1);
			  }
		}
  
}