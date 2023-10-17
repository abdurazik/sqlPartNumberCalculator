using System.Diagnostics;
using CalculateProductNumber;
using Microsoft.VisualBasic;
using NCDK.Charges;



class ConsoleApp
{
    static void Help()
    {
        Console.WriteLine("Help for calcPN:");
        Console.WriteLine("SYNTAX: \"calcPN [arg] [val]\"");
        Console.WriteLine("args:");
        Console.WriteLine("-s\tComputes partnumber for single SMILES string. (eg. -s \"CCOC\")");
        Console.WriteLine("-f\tExpects path to a file where each SMILES string occupies its own line. Output is tab seperated \"SMILES(tab)PartNumber\"");
        Console.WriteLine("NOTE: if no arg is applied, this is the default applied routine");
        Console.WriteLine("-b\tCalls benchmark routine, input being a csv. (eg. -b \"test.csv\")");
        Console.WriteLine("NOTE: csv must have the structure \"PN,SMILES\" WITHOUT headers. Output is tab seperated \"SMILES(tab)PN\"");
        Console.WriteLine("");
        Console.WriteLine("NOTE: outputs as tab-delimated");
        Console.WriteLine("Created by Abdurazik Abdurazik employee of SynquestLaboratories 2023");
    }
    static public async Task BenchMark(string fp)
    {
        Directory.CreateDirectory("Benchmark");
        using (StreamWriter matches = new StreamWriter("Benchmark/match.txt"))
        using (StreamWriter difs = new StreamWriter("Benchmark/dif.txt"))
        using (StreamWriter errors = new StreamWriter("Benchmark/error.txt"))
        {
            
            List<Task> matchTasks = new List<Task>();
            List<Task> difsTasks = new List<Task>();
            matches.Write("PartNumber,CalculatedPartNumber\n");
            difs.WriteLine("PartNumber,CalculatedPartNumber,Smiles\n");
            int matchesBufferBytes= 39 * sizeof(Char);
            int difsBufferBytes= 32 * sizeof(Char);
            
            Stopwatch sw = new Stopwatch();
            StreamReader reader = new StreamReader(fp);
            
            string line;
            sw.Start();
            
            int i = 0; //total count
            int matchCount = 0;
            int difCount = 0;
            int errorCount = 0;
            while ((line = reader.ReadLine()) != null)
            {
                if(difsBufferBytes >= 3500)
                {
                    Task.WhenAll(difsTasks).Wait();
                    difs.Flush();
                }
                if(matchesBufferBytes >= 3500)
                {
                    Task.WhenAll(matchTasks).Wait();
                    matches.Flush();
                }
                
                string[] _s = line.Split(",");
                string pn = _s[0];
                string smiles = _s[1];
                string calculatedPN;
                try
                {
                    calculatedPN = Program.CalculatePN(smiles);
                    if (!pn.Contains(calculatedPN)) 
                    {
                        difsTasks.Add(difs.WriteLineAsync(string.Format("{0}\t{1}\t{2}",pn,calculatedPN,smiles)));
                        difCount ++;
                    }
                    else
                    {
                        matchCount ++;
                        matchTasks.Add(matches.WriteLineAsync(string.Format("{0}\t{1}",pn,calculatedPN)));
                        matchesBufferBytes += (pn.Length + 2 + calculatedPN.Length) * sizeof(Char); 
                    }
                }
                catch (Exception ex)
                {
                    errorCount ++;
                    string errorMsg = ex.Message.Replace("\n",";");
                    
                    errors.WriteLine(string.Format("{0}\t{1}",pn,errorMsg));
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
    static void CalculatePartNumbers(string fp)
    {
        using (StreamWriter @out = new StreamWriter("partnumbers.txt"))
        using (StreamWriter errors = new StreamWriter("errors.txt"))
        {
            StreamReader reader = new StreamReader(fp);
            string? smiles;
            int line =1;
            while ((smiles = reader.ReadLine()) != null)
            {
                line++;
                try
                {
                    string pn = Program.CalculatePN(smiles);
                    @out.WriteLine(string.Format("{0}\t{1}",smiles,pn));
                }
                catch (Exception e)
                {
                    errors.WriteLine("[ERROR] LineNumber:{0} Message: {1}",line,e.Message);
                }
            }
            reader.Close();
        }
    }
    public static async Task Main(string[] args)
    {
        if (args.Contains("-s")) 
        {
            try {Console.WriteLine(Program.CalculatePN(args[1]));}
            catch 
            {
                Console.WriteLine("[ERROR] Please check smiles string, could not parse");
                Console.WriteLine("");
                Help();
            } 
            
        }
        else if(args.Contains("-b")) await BenchMark(args[1]);
        else if(args.Contains("-f")) CalculatePartNumbers(args[1]);
        else
        {
            if (args.Length > 0)
            {
                
            try {CalculatePartNumbers(args[0]);}
            catch 
            {
                Console.WriteLine("[ERROR] : Please check file, Smiles strings must be delimited by new line.");
                Console.WriteLine("###");
                Help();
            }
            }
            else Help();
        }
    }
}