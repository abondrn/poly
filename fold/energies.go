package fold

import (
	"github.com/TimothyStiles/poly/fold/linear/energy_params"
)


const (
	// When calculating the paired minimum free energy if the basepair is isolated,
	// and the seq large, penalize at 1,600 kcal/mol heuristic for speeding this up
	// from https://www.ncbi.nlm.nih.gov/pubmed/10329189
	isolatedBasePairPenalty = 1600
	// maxLenPreCalulated is the length limit of the pre-calculated energies
	// for structures, outside of that range the Jacobson-Stockmayer formula
	// is used see the jacobsonStockmayer() function.
	maxLenPreCalulated = 30

	// minLenForStruct is the minimum length of a nucletide sequence that can
	// create a structure, pentanucleotide sequences form no stable structure
	minLenForStruct = 4

	// loopsAsymmetryPenalty is an energy penalty added for interior loops if
	// the left loop size differs from the right loop size
	// Formula 12 from SantaLucia, 2004
	loopsAsymmetryPenalty = 0.3
	// closingATPenalty is the energy penalty applied to strucfures closed by
	// an AT basepair.
	// Formula 8 from SantaLucia, 2004
	closingATPenalty = 0.5
)

/*
multibranchEnergies holds the a, b, c, d in a linear multi-branch energy

change function.

A - number of helices in the loop
B - number of unpaired nucleotides in the loop
C - coxial stacking in the loop
D - terminal mismatch contributions
E - base composition of the unpaired nucleotides (probably negligible?)
Inferred from:
Supplemental Material: Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Thermodynamics of DNA Structural Motifs, SantaLucia and Hicks, 2004
*/
type multibranchEnergies struct {
	helicesCount, unpairedCount, coaxialStackCount, terminalMismatchCount float64
}

// energy holds two energies, enthaply and entropy
// SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
type energy struct {
	// enthalpy
	enthalpyH float64
	// entropy
	entropyS float64
}

// matchingBasepairEnergy is the energy of matching base pairs
type matchingBasepairEnergy map[string]energy

// loopEnergy is a map[int]Energy where the int is the length of the loop
type loopEnergy map[int]energy

// complementFunc is function to translate a base in its complement
type complementFunc func(rune) rune

// energies holds the needed energy maps, BpEnergy and loopEnergy, to compute
// the folding, it also holds the complement map for the kind of sequence, rna
// or DNA
type energies struct {
	bulgeLoops         loopEnergy
	complement         complementFunc
	danglingEnds       matchingBasepairEnergy
	hairpinLoops       loopEnergy
	multibranch        multibranchEnergies
	internalLoops      loopEnergy
	internalMismatches matchingBasepairEnergy
	nearestNeighbors   matchingBasepairEnergy
	terminalMismatches matchingBasepairEnergy
	triTetraLoops      matchingBasepairEnergy
}


func (e energies) ToEnergyParams(temp float64, scale float64) energy_params.EnergyParams {
	ep := energy_params.EnergyParams{}

	dG := func (nrg energy) int {
		return int(scale*deltaG(nrg.enthalpyH, nrg.entropyS, temp))
	}

	ep.HairpinLoop = make([]int, maxLenPreCalulated)
	for i, nrg := range e.hairpinLoops {
		ep.HairpinLoop[i] = dG(nrg)
	}
	ep.Bulge = make([]int, maxLenPreCalulated)
	for i, nrg := range e.bulgeLoops {
		ep.Bulge[i] = dG(nrg)
	}
	ep.InteriorLoop = make([]int, maxLenPreCalulated)
	for i, nrg := range e.internalLoops {
		ep.InteriorLoop[i] = dG(nrg)
	}

	/* TODO
	// The matrix of free energies for stacked pairs, indexed by the two encoded closing
	// pairs. The list should be formatted as symmetric a `7*7`
	// (`NbDistinguishableBasePairs`*`NbDistinguishableBasePairs`) matrix,
	// conforming to the order explained above. As an example the stacked pair
	// ```
	// 					5'-GU-3'
	// 					3'-CA-5'
	// ```
	// corresponds to the entry StackingPair[1][4] (GC=1, AU=4) which should be
	// identical to StackingPair[4][1] (AU=4, GC=1).
	// size: [NbDistinguishableBasePairs][NbDistinguishableBasePairs]int
	ep.StackingPair

	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	ep.MismatchInteriorLoop [][][]int <- e.internalMismatches
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	ep.Mismatch1xnInteriorLoop [][][]int <- e.internalMismatches
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	ep.Mismatch2x3InteriorLoop [][][]int <- e.internalMismatches
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	ep.MismatchExteriorLoop [][][]int
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	ep.MismatchHairpinLoop [][][]int
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	ep.MismatchMultiLoop [][][]int

	// Energies for the interaction of an unpaired base on the 5' side and
	// adjacent to a helix in multiloops and free ends (the equivalent of mismatch
	// energies in interior and hairpin loops). The array is indexed by the type
	// of pair closing the helix and the unpaired base and, therefore, forms a `7*5`
	// matrix. For example the dangling base in
	// ```
	// 						5'-GC-3'
	// 						3'- G-5'
	// ```
	// corresponds to entry DanglingEndsFivePrime[0][3] (CG=0, G=3).
	//
	// More information about dangling ends: https://rna.urmc.rochester.edu/NNDB/turner04/de.html
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1]int
	ep.DanglingEndsFivePrime [][]int <- e.danglingEnds

	// Same as above for bases on the 3' side of a helix.
	// ```
	// 			       5'- A-3'
	// 			       3'-AU-5'
	// ```
	// corresponds to entry DanglingEndsThreePrime[4][1] (AU=4, A=1).
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1]int
	ep.DanglingEndsThreePrime [][]int <- e.danglingEnds

	// Free energies for symmetric size 2 interior loops. `7*7*5*5` entries.
	// Example:
	// ```
	// 						5'-CUU-3'
	// 						3'-GCA-5'
	// ```
	// corresponds to entry Interior1x1Loop[0][4][4][2] (CG=0, AU=4, U=4, C=2),
	// which should be identical to Interior1x1Loop[4][0][2][4] (AU=4, CG=0, C=2, U=4).
	// size: [NbDistinguishableBasePairs][NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	ep.Interior1x1Loop [][][][]int <- e.interiorLoops

	// Free energies for 2x1 interior loops, where 2 is the number of unpaired
	// nucleotides on the larger 'side' of the interior loop and 1 is the number of
	// unpaired nucleotides on the smaller 'side' of the interior loop.
	// `7*7*5*5*5` entries.
	// Example:
	// ```
	// 						5'-CUUU-3'
	// 						3'-GC A-5'
	// ```
	// corresponds to entry Interior2x1Loop[0][4][4][4][2] (CG=0, AU=4, U=4, U=4, C=2).
	// Note that this matrix is always accessed in the 5' to 3' direction with the
	// larger number of unpaired nucleotides first.
	// size: [NbDistinguishableBasePairs][NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	ep.Interior2x1Loop [][][][][]int <- e.interiorLoops

	// Free energies for symmetric size 4 interior loops. `7*7*5*5*5*5` entries.
	// Example:
	// ```
	// 						5'-CUAU-3'
	// 						3'-GCCA-5'
	// ```
	// corresponds to entry Interior2x2Loop[0][4][4][1][2][2] (CG=0, AU=4, U=4, A=1, C=2, C=2),
	// which should be identical to Interior2x2Loop[4][0][2][2][1][4] (AU=4, CG=0, C=2, C=2, A=1, U=4).
	// size: [NbDistinguishableBasePairs][NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	ep.Interior2x2Loop [][][][][][]int <- e.interiorLoops

	ep.Ninio int
	// size: [NbDistinguishableBasePairs]int
	ep.MultiLoopIntern []int
	ep.MaxNinio int
	*/

	// LogExtrapolationConstant is used to scale energy parameters for hairpin,
	// bulge and interior loops when the length of the loop is greather than
	// `MaxLenLoop`.
	ep.LogExtrapolationConstant = jacobsonStockmayer(1, 1, 0, 1)
	ep.MultiLoopUnpairedNucleotideBonus = 0
	ep.MultiLoopClosingPenalty = 0
	ep.TerminalAUPenalty = 0

	AssignLoop := func(m *map[string]int, k int, loop string, bonus int) {
		if len(loop) == k {
			if *m == nil {
				*m = make(map[string]int)
			}
			(*m)[loop] = bonus
		}
	}

	for loop, nrg := range e.triTetraLoops {
		bonus := dG(nrg)
		AssignLoop(&ep.TriLoop,   5, loop, bonus)
		AssignLoop(&ep.TetraLoop, 6, loop, bonus)
		AssignLoop(&ep.HexaLoop,  8, loop, bonus)
	}

	return ep
}
