using System.Diagnostics;
using CalculateProductNumber;
using Microsoft.VisualBasic;



class ConsoleApp
{
    static public async Task BenchMark(string fp)
    {
        
        using (StreamWriter matches = new StreamWriter("match.csv"))
        using (StreamWriter difs = new StreamWriter("dif.csv"))
        using (StreamWriter errors = new StreamWriter("error.csv"))
        {
            
            List<Task> matchTasks = new List<Task>();
            List<Task> difsTasks = new List<Task>();
            matches.Write("PartNumber,CalculatedPartNumber\n");
            difs.WriteLine("PartNumber,CalculatedPartNumber,Smiles\n");
            int matchesBufferBytes= 39 * sizeof(Char);
            int difsBufferBytes= 32 * sizeof(Char);
            
            Stopwatch sw = new Stopwatch();
            StreamReader reader = new StreamReader(fp);
            reader.ReadLine(); // skip headers
            string line;
            sw.Start();
            
            int i = 0; //total count
            int matchCount = 0;
            int difCount = 0;
            int errorCount = 0;
            while ((line = reader.ReadLine()) != null)
            {
                
                string[] _s = line.Split(",");
                string pn = _s[0];
                string smiles = _s[1];
                string calculatedPN;
                try
                {
                    calculatedPN = Program.CalculatePN(smiles);
                    if (!pn.Contains(calculatedPN)) 
                    {
                        difsTasks.Add(difs.WriteLineAsync(string.Format("{0},{1},{2}",pn,calculatedPN,smiles)));
                        difCount ++;
                    }
                    else
                    {
                        matchCount ++;
                        matchTasks.Add(matches.WriteAsync(string.Format("{0},{1}\n",pn,calculatedPN)));
                        matchesBufferBytes += (pn.Length + 2 + calculatedPN.Length) * sizeof(Char); 
                    }
                }
                catch (Exception ex)
                {
                    errorCount ++;
                    string errorMsg = ex.Message.Replace("\n",";");
                    
                    errors.WriteLineAsync(string.Format("{0},{1}\n",pn,errorMsg));
                    matchesBufferBytes += (pn.Length + 2 + errorMsg.Length) * sizeof(char);
                }
                if (i++ % 20000 == 0) Console.WriteLine("PROCESSED {0} MOLECULES IN {1} SECONDS",i,sw.Elapsed.TotalSeconds);
            }
            sw.Stop();
            Console.WriteLine("\nEND:");
            Console.WriteLine("{0} Errors", errorCount);
            Console.WriteLine("{0} Matching PNs",matchCount);
            Console.WriteLine("{0} Dif PNs",difCount);
            Console.WriteLine("{0} Sturctures processed in {1} seconds processing at {2} molecules/s",i,sw.Elapsed.TotalSeconds,i/sw.Elapsed.TotalSeconds);
        }
    }
    public static async Task Main(string[] args)
    {
        if (args.Contains("-s")) Console.WriteLine(Program.CalculatePN(args[1]));
        else if(args.Contains("-b")) await BenchMark(args[1]);
        
    }
}