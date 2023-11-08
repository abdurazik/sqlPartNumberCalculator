namespace CalculateProductNumber;

using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Text;
using System.Text.RegularExpressions;
using MathNet.Numerics;
using NCDK;
using NCDK.Aromaticities;
using NCDK.Formula;
using NCDK.Graphs;
using NCDK.IO.Formats;
using NCDK.Isomorphisms;
using NCDK.LibIO.CML;
using NCDK.RingSearches;
using NCDK.SMARTS;
using NCDK.Smiles;

public partial class SynquestPartnumberCalculator
{
    static public char FirstCharachter(string smiles,string formula,IAtomContainer molecule)
    {
        bool containsCarbon = false;
        bool containsOxygen = false;
        bool containsNitrogen = false;
        bool Other = false; // Sulfur, Boron, Phosphorus

        if (SmartsPattern.Create("[C,c,N,n,O,o][!#5;!#6;!#7;!#8;!#9;!#17;!#35;!#53;!#14;!#15;!#16;!#17;!#1]").Matches(molecule)) return '9';
        
        
        if(new Regex(@"O(?![a-z])").IsMatch(formula) ) containsOxygen = true;
        if(new Regex(@"N(?![a-z])").IsMatch(formula) ) containsNitrogen = true;
        
        if(
           smiles.Contains("Si")|
           new Regex(@"S(?![a-z])").IsMatch(formula) |
           new Regex(@"P(?![a-z])").IsMatch(formula) |
           new Regex(@"B(?![a-z])").IsMatch(formula)   
           //for smarts aromatic sulfur
           ) Other = true;
        if 
        (
            new Regex(@"C(?![a-z])").IsMatch(formula) 
        ) containsCarbon = true;
        
        if(containsOxygen && containsCarbon)
        {
            if(containsNitrogen & Other) return '8';
            else if (Other) return '6';
            else if (containsNitrogen) return '4';
            else return '2';
        }
        else if(containsNitrogen && containsCarbon)
        {
            if (Other) return '7';
            else return '3';
        }
        else if(Other && containsCarbon) return '5';
        else if (containsCarbon) 
        {
            return '1';
        }
        else return 'M';
    }
    static public char SecondCharachter(string smiles,IAtomContainer molecule)
    {   
        
        if 
        (
            SmartsPattern.Create("*@*@*@*").Matches(molecule) |
            SmartsPattern.Create("C1CC1").Matches(molecule)
        ) //CYCLIC
        {
            if
            (
                SmartsPattern.Create("*@[!#6]@*@*").Matches(molecule) 
            ) return 'H';
            
            else if(new Aromaticity(ElectronDonation.PiBondsModel, Cycles.EdgeShort).Apply(molecule))
            {
                if (SmartsPattern.Create("C=C").Matches(molecule)) //alkene
                {
                    //check if hydrocarbon PN=9 aryl alkene
                    foreach (var substruc in SmartsPattern.Create("C=C").MatchAll(molecule).GetUniqueBonds().ToSubstructures())
                    {
                        IAtom Carbon1 = substruc.Bonds[0].Atoms[0];
                        IAtom Carbon2 = substruc.Bonds[0].Atoms[1];
                        if(HydrocarbonChecker(Carbon1,Carbon2,molecule) && 
                           HydrocarbonChecker(Carbon2,Carbon1,molecule)) return '9';
                    }
                    
                }
                else if (SmartsPattern.Create("C#C").Matches(molecule)) //alkyne
                {
                    //check if hydrocarbon PN=9 aryl alkyne
                    foreach (var substruc in SmartsPattern.Create("C#C").MatchAll(molecule).GetUniqueBonds().ToSubstructures())
                    {
                        IAtom Carbon1 = substruc.Bonds[0].Atoms[0];
                        IAtom Carbon2 = substruc.Bonds[0].Atoms[1];
                        if(HydrocarbonChecker(Carbon1,Carbon2,molecule) && 
                           HydrocarbonChecker(Carbon2,Carbon1,molecule)) return '9';
                    }
                } 
                else if (SmartsPattern.Create("[cR][#6]").Matches(molecule))//alkane
                {
                    if(SmartsPattern.Create("cC(F)(F)F").Matches(molecule)) return '8';
                    //Console.WriteLine("number of substructure {0}",SmartsPattern.Create("c[C;!c][C,H,h]").MatchAll(molecule).GetUniqueBonds().ToSubstructures().Count());
                    foreach (var substruc in SmartsPattern.Create("cC").MatchAll(molecule).GetUniqueBonds().ToSubstructures())
                    {
                        IAtom Carbon1 = substruc.Bonds[0].Atoms[0];
                        IAtom Carbon2 = substruc.Bonds[0].Atoms[1];
                        
                        if
                        (
                            HydrocarbonChecker(Carbon2,Carbon1,molecule) &&
                            HydrocarbonChecker(Carbon1,Carbon2,molecule)
                        ) return '7';
                        
                    }
                    
                }
                return '6';
            }
            else // this is a cycloalkane but need to check for alkenes/ynes as these are dominant structures
            {
                if (Pattern.CreateSubstructureFinder(Chem.MolFromSmiles("C#C")).Matches(molecule)) return '5';
                else if (Pattern.CreateSubstructureFinder(Chem.MolFromSmiles("C=C")).Matches(molecule)) return '4';
                return '2';
            } 
        }
        else //aliphatics
        {
            if (Pattern.CreateSubstructureFinder(Chem.MolFromSmiles("C#C")).Matches(molecule)) return '5';
            else if (Pattern.CreateSubstructureFinder(Chem.MolFromSmiles("C=C")).Matches(molecule)) return '3';
            else return '1';
        }

    }
    static bool HydrocarbonChecker(IAtom atom1, IAtom atom2, IAtomContainer molecule)
    {
        /* This function checks whether the matched aryl alkane/alkene/alkyne is apart a hydrocarbon FG
         used second charachter structure elucidation for Arylalkenes/ynes / Aryl alkanes.
         NOTE this only traverses the half of atom1, reuse function to traverse other half
         of the alkane or alkene or alkyne */
        
        Queue<IAtom> path = new Queue<IAtom>();
        var atoms = from bond in molecule.GetConnectedBonds(atom1)
                    from atom in bond.Atoms
                    where !(atom1.Equals(atom)) && !(atom2.Equals(atom)) && !atom.IsAromatic
                    select atom;
        List<IAtom> searched = new List<IAtom>();
        searched.Add(atom1);
        foreach (IAtom atom in atoms) path.Enqueue(atom);
        //Console.WriteLine("C1 :{0}",atom1);
        //Console.WriteLine("Queued molecules {0}",path.Count);
        while (path.Count > 0)
        {
            IAtom queuedAtom = path.Dequeue();
            searched.Add(queuedAtom);
            if (!new Regex(@"C(?![a-z])|F|C(?=l)|B(?=r)|I|H").IsMatch(queuedAtom.Symbol)) return false; //fixes some bug; lazy fix
            
            //Console.WriteLine("Queued atom: {0}",queuedAtom);
            /* this gets all atoms in bonds not the currently queued atom this will need
            to be further filted for atoms already searched */
            var _atoms = from bond in molecule.GetConnectedBonds(queuedAtom)
                         from __atom in bond.Atoms
                         where !(queuedAtom == __atom)
                         select __atom;
                                    
            
    
            foreach (var _atom in _atoms)
            {
                
                if(_atom.IsAromatic) continue;
                //Console.WriteLine("{0}\nElement:{1} ; Aromatic?:{2}",_atom,_atom.Symbol,_atom.IsAromatic);
                if (!new Regex(@"C(?![a-z])|F|C(?=l)|B(?=r)|I|H").IsMatch(_atom.Symbol)) return false; //Check if atom bonded to polyvalents (other than Carbon)
                
                bool Enqueue = true;
                foreach(var atom in searched) // Check if atom has already been searched
                {
                    if (atom.Equals(_atom))
                    {
                        Enqueue = false;
                        break;
                    }
                }
                if(Enqueue) path.Enqueue(_atom);
                
                
                
            }

        }
        
        return true;

    }
    static public string CnO_FGEvaluator(string smiles,IAtomContainer molecule)
    {
        if
        (
            SmartsPattern.Create("[C,c][C,c](=O)C(=O)O").Matches(molecule)   | //alpha keto-acid/ester
            SmartsPattern.Create("[C,c][C,c](=O)CC(=O)O").Matches(molecule)  | //beta keto-acid/ester
            SmartsPattern.Create("[C,c][C,c](=O)CCC(=O)O").Matches(molecule) | //gamma keto-acid/ester
            SmartsPattern.Create("[C,c]([OH])C(=O)O").Matches(molecule) | //alpha hydroxy-acid/ester
            SmartsPattern.Create("[C,c]([OH])CC(=O)O").Matches(molecule)| //beta hydroxy-acid/ester
            SmartsPattern.Create("[C,c]([OH])CCC(=O)O").Matches(molecule) //gamma hydroxy-acid/ester
        ) return "29";

        else if (SmartsPattern.Create("[C,c](=O)(OC)O[C,c]").Matches(molecule)) return "28"; //carbonates
        else if (SmartsPattern.Create("C(=O)OO").Matches(molecule)) return "27"; //peroxyacids/esters
        else if (SmartsPattern.Create("C(=O)OC(=O)").Matches(molecule)) return "26"; //anhydrides
        else if (SmartsPattern.Create("[O,o]@[C,c](=O)").Matches(molecule)) return "25"; //Lactones TODO:check only if cyclic
        else if 
        (
            SmartsPattern.Create("C=CC(=O)OC").Matches(molecule)  | //alpha unsat-ester
            SmartsPattern.Create("C=CCC(=O)OC").Matches(molecule)   //beta unsat-ester
        )return "24"; //TODO:Check double bond on ester side
        
        else if(SmartsPattern.Create("C(=O)O[C,c]").Matches(molecule)) return "23"; //esters
        else if (SmartsPattern.Create("C(=O)[O-]").Matches(molecule)) return "22"; //carboxylates
        else if(SmartsPattern.Create("C(=O)[OH]").Matches(molecule)) return "21"; //carboxylic acids
        else if
        (
            SmartsPattern.Create("C(=O)C[OH]").Matches(molecule)  |  //alpha hydroxy-ketone
            SmartsPattern.Create("C(=O)CC[OH]").Matches(molecule) |  //beta hydroxy-ketone
            SmartsPattern.Create("C(=O)CCC[OH]").Matches(molecule)   //gamma hydroxy-ketone TODO:check if this should be removed
        ) return "20"; //hydroxy ketones
        
        else if (SmartsPattern.Create("[C,c]=C=O").Matches(molecule)) return "19"; //ketenes
        else if(SmartsPattern.Create("[F,Cl,Br,I][C,c]=O").Matches(molecule)) return "16"; //acid halides
        else if (SmartsPattern.Create("[C,c][C,c](=O)[C,c]").Matches(molecule))
        {
            int uniqueKetones = SmartsPattern.Create("[C,c]C(=O)[C,c]").
                                                MatchAll(molecule).
                                                GetUniqueAtoms().
                                                ToSubstructures().
                                                Count();
            //Console.WriteLine(uniqueKetones);
            if (uniqueKetones > 1) return "18"; //diketones/polycarbonyl
            else return "17"; //ketone
        }
        else if(SmartsPattern.Create("C([OH])([OH])[c,c]").Matches(molecule))return "17" ; //ketone-hydrate (geminal-diol)) return "17";
        else if
        (
            SmartsPattern.Create("[CH]=O").Matches(molecule) |
            SmartsPattern.Create("[CH2]=O").Matches(molecule)
        ) return "15"; //aldehydes
        else if
        (
            SmartsPattern.Create("O-O").Matches(molecule) |
            SmartsPattern.Create("*[O+](*)*").Matches(molecule)
        ) return "14"; //peroxides/oxonium
        
        else if (SmartsPattern.Create("[O-][F,Cl,Br,I]").Matches(molecule)) return "13"; //hypophalite

        int OxygenCount = (from character in smiles
                          where character == 'O'
                          select character).Count();
        if(OxygenCount > 3)
        {
           // Handle potential carbohydrate/crown ethers
        }
        
        if(SmartsPattern.Create("C1OC1").Matches(molecule)) return "09"; //epoxides
        else if(SmartsPattern.Create("[C,c]@[O,o]@[C,c]@[C,c]").Matches(molecule)) return CnO_HetrocycleEval(smiles,molecule); //cyclic ethers
        else if (SmartsPattern.Create("C(OC)(OC)[C,H]").Matches(molecule)) return "08"; //hemiacetal poly ethers handled with Ethers in next evaluation
        else if(SmartsPattern.Create("[c,C][O,o][C,c]").Matches(molecule))
        {
            int etherSubstrucCount = SmartsPattern.Create("[c,C]O[C,c]")
                         .MatchAll(molecule)
                         .GetUniqueAtoms()
                         .ToSubstructures()
                         .Count();
            if(etherSubstrucCount > 1) return "08";
            else return "07";         
        }
        else if(SmartsPattern.Create("[C,c,n]!@O").Matches(molecule)) //alchohols
        {
            int alchoholSubsctructureCount = SmartsPattern.Create("[C,c]O")
                         .MatchAll(molecule)
                         .GetUniqueAtoms()
                         .ToSubstructures()
                         .Count();
            if(alchoholSubsctructureCount > 2) return "05"; //polyols
            else if(alchoholSubsctructureCount ==2) return "04"; //diols

            //If only one alchohol present
            else if(SmartsPattern.Create("C(O)(C)(C)C").Matches(molecule)) return"03"; //tert
            else if(SmartsPattern.Create("[C,c]C(O)[C,c]").Matches(molecule)) return"02"; //second
            else if(SmartsPattern.Create("[c,n]O").Matches(molecule)) return "01"; //weird rule, alchohols bonded to aromatic ring are regarded as primary. Legacy convention
            else return "01"; //must be primary
        }
        return "00";
    }
    static private string CnO_HetrocycleEval(string smiles,IAtomContainer molecule)
    {
        if (SmartsPattern.Create("C(OC)(OC)[C,H]").Matches(molecule)) return "08"; //hemiacetal poly ethers handled with Ethers in next evaluation
        
        
        int etherSubstrucCount = SmartsPattern.Create("[c,C]O[C,c]")
                .MatchAll(molecule)
                .GetUniqueAtoms()
                .ToSubstructures()
                .Count();
        if(etherSubstrucCount > 1) return "08";
        else if (SmartsPattern.Create("[c,C]!@[O]!@[C,c]").Matches(molecule)) return "07";


        
        else if(SmartsPattern.Create("[C,c]!@O").Matches(molecule)) //alchohols
        {
            int alchoholSubsctructureCount = SmartsPattern.Create("[C,c]!@[OH]")
                         .MatchAll(molecule)
                         .GetUniqueAtoms()
                         .ToSubstructures()
                         .Count();
            if(alchoholSubsctructureCount > 2) return "05"; //polyols
            else if(alchoholSubsctructureCount ==2) return "04"; //diols

            //If only one alchohol present
            else if(SmartsPattern.Create("C([OH])(C)(C)C").Matches(molecule)) return"03"; //tert
            else if(SmartsPattern.Create("CC([OH])C").Matches(molecule)) return"02"; //second
            else if(SmartsPattern.Create("[c,n][OH]").Matches(molecule)) return "01"; //weird rule, alchohols bonded to aromatic ring are regarded as primary. Legacy convention
            else if (SmartsPattern.Create("C[OH]").Matches(molecule)) return "01"; //must be primary
        }
        return "10";
    }   
    static public string CnN_FGEvaluator(string smiles,IAtomContainer molecule)
    {

        if (SmartsPattern.CreateSubstructureFinder(Chem.MolFromSmiles("C1=NC=2N=CNC2C=N1")).Matches(molecule)) return "43"; //purines
        
        else if (SmartsPattern.Create("[N-]=[N+]=N[C,c]").Matches(molecule)) return "42"; //azides
        else if
        (
            SmartsPattern.Create("[C,c]C(=N)[N,n]").Matches(molecule) || //amidines
            SmartsPattern.Create("[N,n]C(=N)[N,n]").Matches(molecule)    //guanidines
        )return "41";
        
        else if(SmartsPattern.Create("[C,c]1N=N1").Matches(molecule)) return "40"; //azo
        else if
        (
            SmartsPattern.Create("N-N").Matches(molecule) ||     //Hydrazine
            SmartsPattern.Create("[C,c](=NN)").Matches(molecule) //Hydrazone
        )return "39";
        
        
        else if
        (
            SmartsPattern.Create("C#N").Matches(molecule) || //nitriles
            SmartsPattern.Create("[N+]#[C-]").Matches(molecule) //isonitrile
        )return "37"; 
        
        else if
        (
            SmartsPattern.Create("[C,c]=CN").Matches(molecule) || //enamines
            SmartsPattern.Create("[C,c]=C=N").Matches(molecule)  //ketenamines
        )return "36"; //enamines/ketenamines
        
        else if (SmartsPattern.Create("C=N").Matches(molecule))return "35"; //imines
        else if(SmartsPattern.Create("[C,c]1N[C,c]1").Matches(molecule))return "34"; //aziridines
        else if(SmartsPattern.Create("[!C][N+]([!C])([!C])[!C]").Matches(molecule))return "33"; //ammonium
        else if(SmartsPattern.Create("[C,c]!@N(!@[C,c])[C,c]").Matches(molecule))return "32"; //tertiary amines
        else if(SmartsPattern.Create("[C,c]!@N!@[C,c]").Matches(molecule))return "31"; //secondary amines
        else if(SmartsPattern.Create("[C,c][NH2]").Matches(molecule))return "30";
        
        return "EE";
    }
    static public string CnNO_FGEvaluator(string smiles,IAtomContainer molecule)
    {
        if (SmartsPattern.Create("O1C([NR,nR])CCC1").Matches(molecule)) return "59"; //Nuceleosides
        else if 
        (
           SmartsPattern.Create("C-N").Matches(molecule) &&
           SmartsPattern.Create("[C,c][C,c](=O)O").Matches(molecule)
        ) return aminoAcidHandler(smiles,molecule); 
        
        else if 
        (
            SmartsPattern.Create("[O,o]C(=O)[N,n]").Matches(molecule) || //carbamate
            SmartsPattern.Create("C(=O)NN").Matches(molecule) || //hydrazide
            SmartsPattern.Create("O!=CN").Matches(molecule) ||
            SmartsPattern.Create("ocn").Matches(molecule)
        ) return "56"; //Carbamates
        
        else if 
        (
            SmartsPattern.Create("N[C,c][c,C][N+](=O)[O-]").Matches(molecule)  ||  //alpha nitro-amine
            SmartsPattern.Create("NCC[N+](=O)[O-]").Matches(molecule) ||  //beta  nitro-amine
            SmartsPattern.Create("NCC[N+](=O)[O-]").Matches(molecule) ||  //gamma nitro-amine
            SmartsPattern.Create("C#N[!O]").Matches(molecule)             //nitrile
        ) return "55";
        
        else if 
        (
            SmartsPattern.Create("[N+](=O)[O-]").Matches(molecule) ||
            SmartsPattern.Create("N(=O)=O").Matches(molecule)
        ) return "54"; //nitrated
        else if 
        (
            SmartsPattern.Create("[N,n][C,c](=O)[N,n]").Matches(molecule) || //ureas
            SmartsPattern.Create("C1=CC(=O)NC(=O)N1").Matches(molecule)   || //uracil 
            Pattern.CreateIdenticalFinder(Chem.MolFromSmiles("C1=CC(=O)NC(=O)N1")).Matches(molecule)        
        ) return "53"; 
        
        else if 
        (
            SmartsPattern.Create("[N,n][N,n]=O ").Matches(molecule) || //nitrosoamine
            SmartsPattern.Create("[N,n]=[N+][O-]").Matches(molecule)   //diazene oxide
        ) return "52";
        
        else if (
            SmartsPattern.Create("[C,c]=[N+][O-]").Matches(molecule) || //nitrone
            SmartsPattern.Create("C#[N+][O-]").Matches(molecule) ||    //nitrile oxide
            SmartsPattern.Create("[N+,n+][O-]").Matches(molecule)       //amine oxide
        ) return "50";
        
        else if 
        (
            SmartsPattern.Create("N=C=O").Matches(molecule)    || //isocyanate
            SmartsPattern.Create("O=C=[N-]").Matches(molecule) || //cyanate res
            SmartsPattern.Create("[O+]#C-[N2-]").Matches(molecule) //cyanate res
        ) return "49";
        
        else if 
        (
            SmartsPattern.Create("[C,c](=[O,o])[N,n]").Matches(molecule)  ||  //amide
            SmartsPattern.Create("N@C(=O)").Matches(molecule) ||  //lactam
            SmartsPattern.Create("C(=O)CN").Matches(molecule) ||  //alpha amino-ketone
            SmartsPattern.Create("C(=O)CCN").Matches(molecule) || //beta amino-ketone
            SmartsPattern.Create("C(=O)CCCN").Matches(molecule)   //gamma amino-ketone
        ) return "48";
        
        else if 
        (
            SmartsPattern.Create("C=NO").Matches(molecule) || //oxime
            SmartsPattern.Create("NO").Matches(molecule)      //hydroxylamine
        ) return "47";
        
        else if 
        (
            SmartsPattern.Create("C([OH])C#N").Matches(molecule) || //cyanohydrin
            SmartsPattern.Create("C(=O)C#N").Matches(molecule)
        ) return "46";
        else if 
        (
            SmartsPattern.Create("N=O").Matches(molecule) ||  //nitroso
            SmartsPattern.Create("N[#8;X1v1+0]").Matches(molecule)  //nitroxide
        ) return "45";
        return "EE";
    }
    static public string aminoAcidHandler(string smiles,IAtomContainer molecule)
    {
        if 
        (
            SmartsPattern.Create("[C,c]-[NH][!#1]").Matches(molecule)||
            SmartsPattern.Create("[C,c]-N([!#1])[!#1]").Matches(molecule)
        ) return "58";
        else if 
        (
            SmartsPattern.Create("[C,c][C,c](=O)O[!#1]").Matches(molecule)||
            SmartsPattern.Create("[!#6][C,c](=O)[OH]").Matches(molecule)
        ) return "58";
        return "57";
    }
    static public string Cn_OthersEvaluator(string smiles, IAtomContainer molecule)
    {
        if (smiles.Contains("Si")) return Cn_OtherSilicon(smiles,molecule);
        else if 
        (
            new Regex(@"P(?![ramubod])").IsMatch(smiles) ||
            new Regex(@"(?<!\[C)p").IsMatch(smiles)
        ) return Cn_OtherPhosphorus(smiles,molecule);

        else if 
        (
            new Regex(@"S(?![rgnbem]|c\])").IsMatch(smiles) ||
            new Regex(@"(?<!\[C)s").IsMatch(smiles)
        ) return Cn_OtherSulfur(smiles,molecule);
        
        else if (new Regex(@"B(?![iahrke])").IsMatch(smiles)) return "60";
        else return "EE";
    }
    static public string Cn_OtherSilicon(string smiles, IAtomContainer molecule)
    {
        if 
        (
            SmartsPattern.Create("C(=O)N[Si]").Matches(molecule)|
            SmartsPattern.Create("C=N[Si]").Matches(molecule)   |
            SmartsPattern.Create("[Si]C=N ").Matches(molecule)  |
            SmartsPattern.Create("[Si]C=NO").Matches(molecule)
        ) return "89";
        
        else if (SmartsPattern.Create("[Si]N[Si]").Matches(molecule)) return "88";
        else if (new Regex(@"N(?![abidph])").IsMatch(smiles) | new Regex(@"(?<!\[C)n").IsMatch(smiles)) return "87";
        else if (SmartsPattern.Create("C(=O)O[Si]").Matches(molecule)) return "86";
        else if (SmartsPattern.Create("[Si]O[Si]").Matches(molecule)) return "85";
        else if (SmartsPattern.Create("[Si][O,o][c,n,o]").Matches(molecule)) return "84";
        else if (SmartsPattern.Create("[Si][O,o][C,c]").Matches(molecule)) return "83";
        else if (SmartsPattern.Create("[Si][OH]").Matches(molecule)) return "82";
        else if (SmartsPattern.Create("[Si][F,Cl,Br,I]").Matches(molecule)) return "81";
        else if (SmartsPattern.Create("*[Si](*)(*)*").Matches(molecule)) return "80";
        return "EE";
    }
    static public string Cn_OtherPhosphorus(string smiles,IAtomContainer molecule)
    {
        if 
        (
            SmartsPattern.Create("*-[P,p](=[N,n])(-*)-*").Matches(molecule) ||
            SmartsPattern.Create("[S,O]=PN").Matches(molecule)
        ) return "79";

        else if (SmartsPattern.Create("c1ccccc1P(c1ccccc1)(c1ccccc1)=C ").Matches(molecule)) return "78";
        else if (SmartsPattern.Create("P(=O)(-[O,o]*)(-[o,O]*)[!#8]").Matches(molecule)) return "77";
        else if (SmartsPattern.Create("[Pv4+]").Matches(molecule)) return "76";
        else if (SmartsPattern.Create("P(=O)(-O)(-O)O").Matches(molecule)) return "75";
        else if 
        (
            SmartsPattern.Create("O=P([H,h])(O)O").Matches(molecule) |
            SmartsPattern.Create("O=[PH](O)O").Matches(molecule)
        ) return "74";
        else if (SmartsPattern.Create("[PX3][F,Cl,Br,I]").Matches(molecule)) return "73";
        else if (SmartsPattern.Create("[PX3]").Matches(molecule)) return "72";
        return "EE";
    }
    static public string Cn_OtherSulfur(string smiles,IAtomContainer molecule)
    {
        if (SmartsPattern.Create("[F,Cl,Br,I]S([F,Cl,Br,I])([F,Cl,Br,I])([F,Cl,Br,I])[F,Cl,Br,I]").Matches(molecule)) return "71";
        else if 
        (
            SmartsPattern.Create("*@[S,s]").Matches(molecule) &
            SmartsPattern.Create("*@[!#6;!#16]").Matches(molecule)
        ) return "70";

        else if 
        (
            SmartsPattern.Create("S(=O)(=[O,N])[N,n]").Matches(molecule) ||
            SmartsPattern.Create("S(=O)(=O)[N,n]S(=O)(=O)").Matches(molecule)
        ) return "69";

        else if 
        (
            SmartsPattern.Create("S=C=N").Matches(molecule) ||
            SmartsPattern.Create("SC#N").Matches(molecule)
        ) return "68";
        
        else if (SmartsPattern.Create("S=C").Matches(molecule)) return "67";
        else if (SmartsPattern.Create("[S,s][S,s]").Matches(molecule)) return "66";
        else if (SmartsPattern.Create("[SX3](=O)").Matches(molecule)) return "65";
        else if (SmartsPattern.Create("[SX4](=O)(=O)([O,o])([C,c,H,O,o])").Matches(molecule)) return "64";
        else if (SmartsPattern.Create("[C,c][SX4](=O)(=O)[C,c]").Matches(molecule)) return "63";
        else if (SmartsPattern.Create("S[F,Cl,Br,I]").Matches(molecule)) return "62";
        else if (SmartsPattern.Create("[SX2H]").Matches(molecule)) return "61";
        else if (SmartsPattern.Create("[C,c][SX2][C,c]").Matches(molecule)) return "66";
        return "EE";
    }
    static public string FifthCharachter(string formula)
    {
        bool Hydrogenated = false;
        bool Fluorinated = false;
        bool Chlorinated = false;
        bool Brominated = false;
        bool Iodated = false;

        if (new Regex("H|h(?![a-z])").IsMatch(formula)) Hydrogenated = true;
        if (new Regex("F(?![a-z])").IsMatch(formula)) Fluorinated = true;
        if (new Regex("Cl").IsMatch(formula)) Chlorinated = true;
        if (new Regex("Br").IsMatch(formula)) Brominated = true;
        if (new Regex("I(?![a-z])").IsMatch(formula)) Iodated = true;

        if (Fluorinated)
        {
            if (Iodated)
            {
                if (Hydrogenated && Chlorinated && Brominated) return "Y-"; 
                else if (Chlorinated && Brominated) return "X-"; 
                else if (Hydrogenated && Brominated) return "U-"; 
                else if (Brominated ) return "T-"; 
                else if (Hydrogenated && Chlorinated ) return "P-"; 
                else if (Chlorinated) return "N-"; 
                else if (Hydrogenated ) return "K-"; 
                else return "J-"; 
            }
            else if (Brominated)
            {
                if (Hydrogenated && Chlorinated && Brominated) return "F-"; 
                else if (Chlorinated && Brominated) return "E-"; 
                else if (Hydrogenated && Brominated) return "B-"; 
                else return "A-"; 
            }
            if (Chlorinated)
            {
                if (Hydrogenated && Chlorinated) return "7-"; 
                else return "6-"; 
            }
            if (Hydrogenated) return "3-"; 
            return "2-"; 
        }
        else if (Chlorinated)
        {
            if (Iodated)
            {
                if(Brominated && Hydrogenated) return "W-";
                if(Brominated) return "V-";
                else if (Hydrogenated) return "M-";
                else return "L-";
            }
            else if (Brominated)
            {
                if (Hydrogenated) return "D-";
                else return "C-";
            } 
            else if (Hydrogenated) return "5-";
            else return "4-";
        }
        else if (Brominated)
        {
            if (Iodated && Hydrogenated) return "S-";
            else if (Iodated) return "R-";
            else if (Hydrogenated) return "9-";
            else return "8-";
        }
        else if (Iodated)
        {
            if (Hydrogenated) return "H-";
            else return "G-";
        }
        else if (Hydrogenated) return "1-";
        else return "0-";
    }
    static public string CalculatePN(string smiles)
    {

        if (smiles == null || smiles == "") return "EEEE-";
        
        char firstCharachter = 'E';
        char secondCharachter ='E';
        string FGCode = "EE";
        int FGCodeint = 0; //used only for PN's that start with 9 or M
        string fifthCharachter = "E-";
        
        
        IAtomContainer molecule = Chem.MolFromSmiles(smiles);
        string inchi = Chem.MolToInChI(molecule);
        string formula = inchi.Split("/")[1];
        firstCharachter = FirstCharachter(smiles,formula,molecule);
        secondCharachter = SecondCharachter(smiles,molecule);
        fifthCharachter = FifthCharachter(formula);
        
        switch(firstCharachter)
        {
            case '1':
                FGCode = "00";
                break;
            
            case '2':
                FGCode = CnO_FGEvaluator(smiles,molecule);
                break;
            
            case '3':
                FGCode = CnN_FGEvaluator(smiles,molecule);
                if (FGCode == "EE") 
                {
                    if 
                    (
                        SmartsPattern.Create("n").Matches(molecule) ||
                        SmartsPattern.Create("C@N(C)@C").Matches(molecule) 
                    ) FGCode =  "32";
                    else FGCode = "31";
                }
                break;
            
            case '4':
                FGCode = CnNO_FGEvaluator(smiles,molecule);
                if (secondCharachter == 'H')
                {
                    if(FGCode == "00" || FGCode == "EE") FGCode = CnN_FGEvaluator(smiles,molecule);        
                    
                    if 
                    (
                        SmartsPattern.Create("[O!R,o!r]").Matches(molecule) &&
                        (FGCode == "00" || FGCode == "EE")
                    )
                    {
                        FGCode = CnO_FGEvaluator(smiles,molecule);
                    }
                    if (FGCode == "00" || FGCode == "EE")
                    {
                         if 
                        (
                            SmartsPattern.Create("n").Matches(molecule) ||
                            SmartsPattern.Create("C@N(C)@C").Matches(molecule) 
                        ) FGCode =  "32";
                        else FGCode = "31";
                    }
                }
                else
                {
                    if (FGCode == "00" || FGCode == "EE") FGCode = CnN_FGEvaluator(smiles,molecule);
                    else if (FGCode == "00" || FGCode == "EE") FGCode = CnO_FGEvaluator(smiles,molecule);
                }
                
                break;
            
            case '5':
            case '6':
            case '7':
            case '8':
                FGCode = Cn_OthersEvaluator(smiles,molecule);
                break;
            
            case '9':
            case 'M':
                if (firstCharachter == 'M') secondCharachter = '0';
                int largestAN = 0;
                int largestHalogen = 0;
                foreach (var atom in molecule.Atoms)
                {
                    if 
                    (
                        atom.AtomicNumber != 6 &&
                        atom.AtomicNumber != 9 &&
                        atom.AtomicNumber != 17 &&
                        atom.AtomicNumber != 35 &&
                        atom.AtomicNumber != 53
                    )
                    {
                        if (atom.AtomicNumber > largestAN) largestAN = atom.AtomicNumber;
                    }

                    if 
                    (
                        atom.AtomicNumber == 9  ||
                        atom.AtomicNumber == 17 ||
                        atom.AtomicNumber == 35 ||
                        atom.AtomicNumber == 53
                    )
                    {
                        if (atom.AtomicNumber > largestHalogen) largestHalogen = atom.AtomicNumber;
                    }
                    
                    

                }
                
                if (largestAN == 1 || largestAN == 0)
                    {
                        if (largestHalogen.ToString().Length == 1) FGCode = string.Format("0{0}",largestHalogen);
                        else FGCode = largestHalogen.ToString();
                    }
                    else
                    {
                        if (largestAN.ToString().Length == 1) FGCode = string.Format("0{0}",largestAN);
                        else FGCode = largestAN.ToString();
                    }
                break;
            
            default:
                FGCode = "EE";
                break;
        }
        
        return string.Format("{0}{1}{2}-{3}",firstCharachter,secondCharachter,FGCode,fifthCharachter);
    }
    

}
