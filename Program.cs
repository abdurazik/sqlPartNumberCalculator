namespace CalculateProductNumber;

using System.ComponentModel.Design;
using System.Net.Mail;
using System.Runtime.InteropServices;
using System.Text.RegularExpressions;
using Microsoft.VisualBasic;
using NCDK;
using NCDK.Aromaticities;
using NCDK.Graphs;
using NCDK.IO;
using NCDK.Isomorphisms;
using NCDK.QSAR.Descriptors.Moleculars;
using NCDK.SMARTS;

class Program
{
    static public char FirstCharachter(string smiles, IAtomContainer molecule)
    {
        

        bool containsOxygen = false;
        bool containsNitrogen = false;
        bool Other = false; // Sulfur, Boron, Phosphorus

        if(new Regex(@"O(?![a-z])").IsMatch(smiles)) containsOxygen = true;
        if(new Regex(@"N(?![a-z])").IsMatch(smiles)) containsNitrogen = true;
        if(new Regex(@"S(?![a-z])").IsMatch(smiles) |
           new Regex(@"P(?![a-z])").IsMatch(smiles) |
           new Regex(@"B(?![a-z])").IsMatch(smiles)) Other = true;
        
        if(containsOxygen)
        {
            if(containsNitrogen & Other) return '8';
            else if (Other) return '6';
            else if (containsNitrogen) return '4';
            else return '2';
        }
        else if(containsNitrogen)
        {
            if (Other) return '7';
            else return '3';
        }
        else if(Other) return '5';
        return '1'; //Carbon only if all are false
    }
    static public char SecondCharachter(string smiles,IAtomContainer molecule)
    {
        bool isCyclic = Cycles.FindAll(molecule).GetNumberOfCycles() > 0;
        
        if (isCyclic)
        {
            if(SmartsPattern.Create("[!CR;!cR]").Matches(molecule)) return 'H';
            else if(new Aromaticity(ElectronDonation.PiBondsModel, Cycles.AllSimpleFinder).Apply(molecule))
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
                else if (SmartsPattern.Create("cC").Matches(molecule))//alkane
                {
                     foreach (var substruc in SmartsPattern.Create("cC").MatchAll(molecule).GetUniqueBonds().ToSubstructures())
                    {
                        IAtom Carbon1 = substruc.Bonds[0].Atoms[0];
                        IAtom Carbon2 = substruc.Bonds[0].Atoms[1];
                        if(HydrocarbonChecker(Carbon1,Carbon2,molecule) && 
                           HydrocarbonChecker(Carbon2,Carbon1,molecule)) return '7';
                    }
                }
                return '6';
            }
            else // this is a cycloalkane but need to check for alkenes/ynes as these are dominant structures
            {
                if (Pattern.CreateSubstructureFinder(Chem.MolFromSmiles("C#C")).Matches(molecule)) return '5';
                else if (Pattern.CreateSubstructureFinder(Chem.MolFromSmiles("C=C")).Matches(molecule)) return '3';
                return '4';
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
                    where !(atom1.Equals(atom)) && !(atom2.Equals(atom))
                    select atom;
        List<IAtom> searched = new List<IAtom>();
        searched.Add(atom1);
        foreach (IAtom atom in atoms) path.Enqueue(atom);
        
        while (path.Count > 0)
        {
            IAtom queuedAtom = path.Dequeue();
            searched.Add(queuedAtom);
            
            
            /* this gets all atoms in bonds not the currently queued atom this will need
            to be further filted for atoms already searched */
            var _atoms = from bond in molecule.GetConnectedBonds(queuedAtom)
                         from __atom in bond.Atoms
                         where !(queuedAtom == __atom)
                         select __atom;
                                    
            
    
            foreach (var _atom in _atoms)
            {
                
                if(_atom.IsAromatic) continue;
                if (!(new Regex(@"C(?![a-z])|F|C(?=l)|B(?=r)|I").IsMatch(_atom.Symbol))) return false; //Check if atom bonded to polyvalents (other than Carbon)
                
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
    static void Main(string[] args)
    {
        string smiles = "C2CCCCC2CCCCC=C";
        
        IAtomContainer molecule = Chem.MolFromSmiles(smiles);
        
        Console.WriteLine(SecondCharachter(smiles,molecule));

        
    }
}
