/*
 * reads files and returns string array
 */


import java.io.*;
import java.util.*;                    // for ArrayList
import java.util.regex.*;              // for regular expression methods (s. also RegexTestHarness.java)


public class ReadFile
{   private File inputFile;

    public ReadFile()
    { 	this(new File(""));
    }

    public ReadFile(File f)
    {  inputFile = f;
    }

    public HashMap<String,Integer> getLengthMap(int size)
		  {  HashMap<String,Integer> slMap = new HashMap<String,Integer>(size);

		     try
		     {  BufferedReader in = new BufferedReader(new FileReader(this.inputFile)); 			      // no column titles
          String currentLine = in.readLine();


          // check if fasta
          Pattern pFasta = Pattern.compile(">");
          boolean fasta=false;
          Matcher mFasta = pFasta.matcher(currentLine);
						    if (mFasta.find())
						      fasta = true;
						    String id="", seq="";
						    if (fasta)
						    {  while (currentLine!=null)
			          {  mFasta = pFasta.matcher(currentLine);
							         if (mFasta.find())
							         {  id = currentLine.substring(1);  // id
							         }
							         else
							         {  seq = currentLine.trim();
							            slMap.put(id,seq.length());
							         }
							         currentLine = in.readLine();
						       }
						    }
						    else
						    {  Pattern pTab = Pattern.compile("\t");
						    	  while (currentLine!=null)
			          { 	Matcher mTab = pTab.matcher(currentLine);
									       if (mTab.find())
									       {  String[] item = pTab.split(currentLine);
									          id = item[0];
									          seq = item[1];
									          slMap.put(id,seq.length());
									       }
									       currentLine = in.readLine();
						       }
						    }

						    in.close();
						    System.out.println("Input file " + inputFile.getName() + " read.");
		     }
				   catch (Exception exc)
	      { System.out.println("Error: " + exc);
	        System.exit(1);
	      }
       return slMap;
				}

    public ArrayList<String[]> getStringArray()   // add delimiter option? (default: tab)
		  {  ArrayList<String[]> sList = new ArrayList<String[]>(10);

		     try
		     {  BufferedReader in = new BufferedReader(new FileReader(this.inputFile)); 			      // no column titles
          String currentLine = in.readLine();
          Pattern pTab = Pattern.compile("\t");
          Pattern pFasta = Pattern.compile(">");
          boolean fasta=false;
          Matcher mFasta = pFasta.matcher(currentLine);
						    if (mFasta.find())
						      fasta = true;
						    if (fasta)
						    {  int counter=0;
             String[] seqRow = new String[2];
             seqRow[1]="";
						    	  while (currentLine!=null)
			          {  mFasta = pFasta.matcher(currentLine);
							         if (mFasta.find())
							         {  if (counter>0)
							            {  sList.add(seqRow);
							               seqRow = new String[2];
							               seqRow[1]="";
							            }
							            seqRow[0] = currentLine.substring(1);  // id
							         }
							         else
							         {  seqRow[1] = seqRow[1] + currentLine.trim();
							         }
							         counter++;
							         currentLine = in.readLine();
						       }
						       sList.add(seqRow);  // last entry
						    }
						    else
						    {  while (currentLine!=null)
			          { 	Matcher mTab = pTab.matcher(currentLine);
									       if (mTab.find())
									       {  String[] item = pTab.split(currentLine);
									          sList.add(item);
									       }
									       else
									       {  String[] currentColumn = new String[1];
									          currentColumn[0] = currentLine.trim();
									          if (!currentColumn[0].equals(""))
									            sList.add(currentColumn); // if only one column (-> no tabs!)
									       }
									       currentLine = in.readLine();
						       }
						    }

						    in.close();
						    System.out.println("Input file " + inputFile.getName() + " read.");
		     }
				   catch (Exception exc)
	      { System.out.println("Error: " + exc);
	        System.exit(1);
	      }
       return sList;
				}

				public ArrayList<String[]> getStringArray(String separator, boolean ignoreFirstRow)   // add delimiter option? (default: tab)
		  {  ArrayList<String[]> sList = new ArrayList<String[]>(10);

		     try
		     {  BufferedReader in = new BufferedReader(new FileReader(this.inputFile)); 			      // no column titles
          String currentLine = in.readLine();
          if (ignoreFirstRow)
          {  currentLine = in.readLine();
          }
          Pattern pSep = Pattern.compile(separator);
          while (currentLine!=null)
			       {  Matcher mSep = pSep.matcher(currentLine);
						       if (mSep.find())
						       {  String[] item = pSep.split(currentLine);//System.out.println(currentLine);System.out.println("item.length " + item.length);
						          String[] itemNew = new String[item.length-1];
						          if (item[0].equals(""))						       // in case item[0] is empty -> remove
						          {  for (int i=1; i<item.length; i++)
						             {  //System.out.println("item"+i+" "+item[i]);
						                itemNew[i-1] = item[i]; //System.out.println("itemNew"+i+" "+itemNew[i-1]);
						             }
						          }
						          else
						            itemNew = item;
						          sList.add(itemNew);
						       }
						       else
						       {  String[] currentColumn = new String[1];
						          currentColumn[0] = currentLine.trim();
						          if (!currentColumn[0].equals(""))
						            sList.add(currentColumn); // if only one column (-> no tabs!)
						       }
						       currentLine = in.readLine();
						    }

						    in.close();
						    System.out.println("Input file " + inputFile.getName() + " read.");
		     }
				   catch (Exception exc)
	      { System.out.println("Error: " + exc);
	        System.exit(1);
	      }
       return sList;
				}

				public ArrayList<String[]> getHeaderArray()   // add delimiter option? (default: tab)
		  {  ArrayList<String[]> sList = new ArrayList<String[]>(10);

		     try
		     {  BufferedReader in = new BufferedReader(new FileReader(this.inputFile)); 			      // no column titles
          String currentLine = in.readLine();
          Pattern pTab = Pattern.compile("\t");
          Matcher mTab = pTab.matcher(currentLine);
			       if (mTab.find())
			       {  String[] item = pTab.split(currentLine);
			          sList.add(item);
			       }
			       else
			       {  String[] currentColumn = new String[1];
			          currentColumn[0] = currentLine.trim();
			          if (!currentColumn[0].equals(""))
			            sList.add(currentColumn); // if only one column (-> no tabs!)
			       }

						    in.close();
		     }
				   catch (Exception exc)
	      { System.out.println("Error: " + exc);
	        System.exit(1);
	      }
       return sList;
				}

				public ArrayList<String[]> getStringArray(boolean withHeader)
		  {  ArrayList<String[]> sList = new ArrayList<String[]>(10);

		     try
		     {  BufferedReader in = new BufferedReader(new FileReader(this.inputFile)); 			      // no column titles
          String currentLine = in.readLine();
          if (withHeader)
          {  currentLine = in.readLine();
          }
          Pattern pSep = Pattern.compile("\t");
          while (currentLine!=null)
			       {  Matcher mSep = pSep.matcher(currentLine);
						       if (mSep.find())
						       {  String[] item = pSep.split(currentLine);//System.out.println(currentLine);System.out.println("item.length " + item.length);
						          String[] itemNew = new String[item.length-1];
						          if (item[0].equals(""))						       // in case item[0] is empty -> remove
						          {  for (int i=1; i<item.length; i++)
						             {  //System.out.println("item"+i+" "+item[i]);
						                itemNew[i-1] = item[i]; //System.out.println("itemNew"+i+" "+itemNew[i-1]);
						             }
						          }
						          else
						            itemNew = item;
						          sList.add(itemNew);
						       }
						       else
						       {  String[] currentColumn = new String[1];
						          currentColumn[0] = currentLine.trim();
						          if (!currentColumn[0].equals(""))
						            sList.add(currentColumn); // if only one column (-> no tabs!)
						       }
						       currentLine = in.readLine();
						    }

						    in.close();
						    System.out.println("Input file " + inputFile.getName() + " read.");
		     }
				   catch (Exception exc)
	      { System.out.println("Error: " + exc);
	        System.exit(1);
	      }
       return sList;
				}


	   public String getFileNameBase()        // for finding filename without extension
	   {  StringBuffer name = new StringBuffer();
	      StringBuffer nameBase = new StringBuffer();
	      name.append(this.inputFile.getName());
	      for (int i=0; i<name.length(); i++)
	      {  String c = name.charAt(i)+"";            // (conversion to String necessary, otherwise "equals"-operation not possible)
	         if (c.equals("."))
	           break;
	         nameBase.append(name.charAt(i));
	      }
	      return nameBase.toString();
	   }

     public String getName()
	   {  return this.inputFile.getName();
	   }

    public String getPath()
	   {  return this.inputFile.getPath();
	   }

	   public String getCanonPath()
	   {  String path = "";
	      try
	      {	 path = this.inputFile.getCanonicalPath();
	      }
	      catch (Exception exc)
	      { System.out.println("Error: " + exc);
	        System.exit(1);
	      }
	      return path;
	   }

	   public String getCanonParent()
	   {  String path = "", parent="";
	      try
	      {	 path = this.inputFile.getCanonicalPath();
	         parent = path.substring(0,path.length()-this.inputFile.getName().length()-1);
	      }
	      catch (Exception exc)
	      { System.out.println("Error: " + exc);
	        System.exit(1);
	      }
	      return parent;
	   }

	   public String getParent()
	   {  return this.inputFile.getParent();
	   }

	   public String getParent(String path)
	   {  return path.substring(0,(int)(path.length()-this.inputFile.getName().length()));
	   }
}