package fold

import (
	"fmt"
	"math"
	"strings"

	"github.com/TimothyStiles/poly/checks"
)


// subsequence represent an interval of bases in the sequence that can contain
// a inward structure.
type subsequence struct {
	start, end int
}

// nucleicAcidStructure is single structure with a free energy, description, and inward children
type nucleicAcidStructure struct {
	// description is the description of the nucleicAcidStructure
	description string
	// inner is list of inner structure represented as intervals in the
	// sequence
	inner []subsequence
	// energy is the energy of the nucleicAcidStructure
	energy float64
}

// Equal returns true if two nucleicAcidStructures are equal
func (structure nucleicAcidStructure) Equal(other nucleicAcidStructure) bool {
	if len(structure.inner) != len(other.inner) {
		return false
	}
	for i, val := range structure.inner {
		if val != other.inner[i] {
			return false
		}
	}
	return structure.energy == other.energy
}

// Valid returns true if the NucleicAcidStructure is valid
func (structure nucleicAcidStructure) Valid() bool {
	return structure.energy != math.Inf(1) && structure.energy != math.Inf(-1)
}

// defaultStructure is the default (zero value) nucleic acid structure, it used
// mostly to initialize the caches, see Context
var defaultStructure = nucleicAcidStructure{
	description: "",
	energy:      math.Inf(-1),
}

// invalidStructure represent an invalid nucleic acid structure
var invalidStructure = nucleicAcidStructure{
	description: "",
	energy:      math.Inf(1),
}

// context holds the energy caches, energy maps, sequence, and temperature
// needed in order to compute the folding energy and structures.
type context struct {
	energies                   energies
	seq                        string
	pairedMinimumFreeEnergyV   [][]nucleicAcidStructure
	unpairedMinimumFreeEnergyW [][]nucleicAcidStructure
	temp                       float64
}

// newFoldingContext returns a context ready to use, in case of error
// the returned FoldingContext is empty.
func newFoldingContext(seq string, temp float64) (context, error) {
	seq = strings.ToUpper(seq)

	// figure out whether it's DNA or rna, choose energy map
	var energyMap energies
	switch {
	case checks.IsDNA(seq):
		energyMap = dnaEnergies
	case checks.IsRNA(seq):
		energyMap = rnaEnergies
	default:
		return context{}, fmt.Errorf("the sequence %s is not RNA or DNA", seq)
	}

	var (
		sequenceLength = len(seq)
		vCache         = make([][]nucleicAcidStructure, sequenceLength)
		wCache         = make([][]nucleicAcidStructure, sequenceLength)
		row            = make([]nucleicAcidStructure, sequenceLength)
	)
	for nucleicAcidIndex := 0; nucleicAcidIndex < sequenceLength; nucleicAcidIndex++ {
		row[nucleicAcidIndex] = defaultStructure
	}
	for j := 0; j < sequenceLength; j++ {
		vCache[j] = make([]nucleicAcidStructure, sequenceLength)
		copy(vCache[j], row)

		wCache[j] = make([]nucleicAcidStructure, sequenceLength)
		copy(wCache[j], row)
	}
	ret := context{
		energies:                   energyMap,
		seq:                        seq,
		pairedMinimumFreeEnergyV:   vCache,
		unpairedMinimumFreeEnergyW: wCache,
		temp:                       temp + 273.15, // kelvin
	}

	// fill the cache
	_, err := unpairedMinimumFreeEnergyW(0, sequenceLength-1, ret)
	if err != nil {
		return context{}, fmt.Errorf("error filling the caches for the FoldingContext: %w", err)
	}
	return ret, nil
}

// Result holds the resulting structures of the folded s
type Result struct {
	structs []nucleicAcidStructure
}

// DotBracket returns the dot-bracket notation of the secondary nucleic acid
// structure resulting from folding a sequence.
//
// Dot-bracket notation, consisting in a balanced parentheses string composed
// by a three-character alphabet {.,(,)}, that can be unambiguously converted
// in the RNA secondary structure. See example_test.go for a small example.
func (r Result) DotBracket() string {
	if len(r.structs) == 0 {
		return ""
	}
	lastStructEnd := 0
	for _, structure := range r.structs {
		for _, innerSubsequence := range structure.inner {
			if innerSubsequence.end > lastStructEnd {
				lastStructEnd = innerSubsequence.end
			}
		}
	}
	lastStructEnd += 1
	result := make([]byte, lastStructEnd)
	for i := range result {
		result[i] = '.'
	}
	for _, structure := range r.structs {
		if len(structure.inner) == 1 {
			innerSubsequence := structure.inner[0]
			result[innerSubsequence.start] = '('
			result[innerSubsequence.end] = ')'
		}
	}
	return string(result)
}

// MinimumFreeEnergy return just the delta G of the structures resulting from
// folding a sequence.
//
// Returns the minimum free energy of the folded sequence
func (r Result) MinimumFreeEnergy() float64 {
	if len(r.structs) == 0 {
		// invalid
		return math.Inf(1)
	}

	summedEnergy := 0.0
	for _, structure := range r.structs {
		summedEnergy += structure.energy
	}
	return summedEnergy
}
