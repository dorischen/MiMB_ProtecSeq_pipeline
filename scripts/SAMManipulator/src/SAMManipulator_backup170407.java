/**
 * @(#)SAMManipulator.java
 *
 * input: SAM alignment file, read ids
 * output: SAM alignment file with or without specified reads
 

 *
 * @doris chen
 * @version 1.00 2016/06/17
 */

import java.io.*;
import java.util.*;                    // for ArrayList
import java.util.regex.*;              // for regular expression methods (s. also RegexTestHarness.java)


public class SAMManipulator
{  static File inputFileS, inputFileR;
   static String option, outfileBase;  

   public static void main (String args[]) throws Exception
   {
      // default
      option="";
      outfileBase="";

      try
      {  if (args.length>=4)
			      {  Pattern pSam = Pattern.compile("^-sam");
			         Pattern pReads = Pattern.compile("^-reads");
			         Pattern pOption = Pattern.compile("^-option");
			         Pattern pOut = Pattern.compile("^-out");
			        
            for (int i=0; i<args.length; i++)
			         {  Matcher mSam = pSam.matcher(args[i]);
			            if (mSam.find())
			            {  inputFileS = new File(args[i+1]);
			               if ( !inputFileS.exists() ) throw new Exception();
			            }
               else
               {  Matcher mReads = pReads.matcher(args[i]);
									            if (mReads.find())
									            {   inputFileR = new File(args[i+1]);
			                      if ( !inputFileR.exists() ) throw new Exception();
									            }
									            else
									            {  Matcher mOption = pOption.matcher(args[i]);
												            if (mOption.find())
												            {  option = args[i+1];
												               if ( !option.equals("delete") && !option.equals("keep") ) throw new Exception();
												            }  else
												            {  Matcher mOut = pOut.matcher(args[i]);
															            if (mOut.find())
															            {  outfileBase = args[i+1];															               
															            }
															         }
												         }
															 }
			         }

            process_sam(inputFileS, inputFileR, option, outfileBase);

			      }
			      else throw new Exception();
			   }
			   catch (Exception exc)
			   {  System.out.println("");
			      System.out.println(">java SAMManipulator -sam SAM_FILE -reads READ_FILE -option delete|keep -out OUTFILE_BASE");
			      System.out.println("ad READ_FILE: read ids, one line each");
			      System.out.println("ad option .. delete -> reads deleted from sam; keep -> reads kept (rest of reads deleted)");
			      System.out.println("ad OUTFILE_BASE: file name for resulting sam file (.sam appended)");
			      	
         System.out.println("Error: " + exc.toString() + '\n');
			      System.out.println("");
			   }

   }
 

	  public static void process_sam(File inputFileS, File inputFileR, String option, String outfileBase)
	  {  try
			   {  boolean delete=false;
			      String suffix = "retained";
			   	  if(option.equals("delete"))
			   	  {  delete = true;
			   	     suffix = "deleted";
			   	  }
			   	  
			   	  // get reads
         ReadFile rfR = new ReadFile(inputFileR);
         BufferedReader inR = new BufferedReader(new FileReader(inputFileR));
         
         ArrayList<String> readList = new ArrayList(10000);
         String currentLineR = inR.readLine();
         while(currentLineR!=null)
         {  readList.add(currentLineR.trim());
         	  currentLineR = inR.readLine();
         }
         inR.close();

         // go through sam file
         ReadFile rfS = new ReadFile(inputFileS);
         Pattern pTab = Pattern.compile("\t"); 
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
         while(currentLineS!=null)
         {  // go through rows, split into items
            String[] row = pTab.split(currentLineS);
            currentRead = row[0];
            if(readList.contains(currentRead))
            {  if(!delete)
            	    outS.println(currentLineS.trim());            	
            }
            else
            {  if(delete)
            	    outS.println(currentLineS.trim()); 
            }
            
            currentLineS = inS.readLine();
			      }
         inS.close();
         outS.close();
         
         System.out.println(readList.size()+" reads "+suffix+".");
         System.out.println("");

     }
			  catch (Exception exc)
			  { System.out.print("Error: " + exc);
			    System.exit(1);
			  }
		}
  
}