package annotate

import (
	"encoding/csv"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"

	"github.com/TimothyStiles/poly/io/fasta"
	"gopkg.in/yaml.v3"
)

// Task represents a bioinformatics task.
type Task interface {
	// Run takes paths to database/input/output files, and logfile log
	Run(dbPath, inPath, outPath string, log *os.File) error
}

// BlastTask represents the BLAST task.
type BlastTask struct{}

func CreateTempFasta(seq string) (string, error) {
	bytes, err := fasta.Build([]fasta.Fasta{{"temp", seq}})
	if err != nil {
		return "", err
	}

	query, err := os.CreateTemp("", "query_*.fasta")
	if err != nil {
		return "", err
	}
	_, err = query.Write(bytes)
	if err != nil {
		return "", err
	}
	query.Close()
	return query.Name(), nil
}

func (t BlastTask) Run(dbPath, inPath, outPath string, log *os.File) error {
	flags := "qstart qend sseqid sframe pident slen qseq length sstart send qlen evalue"
	cmd := exec.Command("blastn", "-task", "blastn-short", "-query", inPath, "-out", outPath,
		"-db", dbPath, "-outfmt", fmt.Sprintf("6 %s", flags))
	cmd.Stdout = log
	cmd.Stderr = log

	return cmd.Run()
}

// DiamondTask represents the Diamond task.
type DiamondTask struct{}

func (t DiamondTask) Run(dbPath, inPath, outPath string, log *os.File) error {
	flags := "qstart qend sseqid pident slen qseq length sstart send qlen evalue"
	cmd := exec.Command("diamond", "blastx", "-d", dbPath, "-q", inPath, "-o", outPath,
		"--outfmt", fmt.Sprintf("6 %s", flags))
	cmd.Stdout = log
	cmd.Stderr = log

	return cmd.Run()
}

// InfernalTask represents the Infernal task.
type InfernalTask struct{}

func (t InfernalTask) Run(dbPath, inPath, outPath string, log *os.File) error {
	flags := "--cut_ga --rfam --noali --nohmmonly --fmt 2"
	cmd := exec.Command("cmscan", flags, "--tblout", outPath, "--clanin", dbPath, inPath)
	cmd.Stdout = log
	cmd.Stderr = log

	return cmd.Run()
}

func parseInfernal(filename, seq string) error {
	// Your Infernal parsing logic goes here.
	return nil
}

type Database struct {
	Version    string   `yaml:"version"`
	Method     string   `yaml:"method"`
	Location   string   `yaml:"location"`
	Priority   int      `yaml:"priority"`
	Parameters []string `yaml:"parameters"`
	Details    struct {
		DefaultType string `yaml:"default_type"`
		Location    string `yaml:"location"`
		Compressed  bool   `yaml:"compressed"`
	} `yaml:"details"`
}

func LoadDatabases(path string) (Databases, error) {
	unparsed, err := os.Open(path)
	if err != nil {
		return Databases{}, err
	}

	// Define a map to store the unmarshaled data
	var parsed Databases

	// Unmarshal YAML data into the map
	err = yaml.NewDecoder(unparsed).Decode(&parsed)
	if err != nil {
		return Databases{}, err
	}

	dir := filepath.Dir(path)
	for _, v := range parsed {
		v.Location = filepath.Join(dir, v.Location)
	}

	return parsed, nil
}

type Hit map[string]string

type Databases map[string]Database

func getRawHits(query string, linear bool, dbs Databases) ([]Hit, error) {
	logFile, err := os.Create("output.log")
	if err != nil {
		return []Hit{}, err
	}
	defer logFile.Close()

	inPath, err := CreateTempFasta(query)
	if err != nil {
		return []Hit{}, err
	}
	all_hits := []Hit{}
	for name, db := range dbs {
		hits, err := Blast(inPath, name, db, logFile)
		if err != nil {
			return []Hit{}, err
		}

		all_hits = append(all_hits, hits...)
	}
	err = os.Remove(inPath)
	return []Hit{}, err
}

func readCSV(filename string) ([]Hit, error) {
	// Open the CSV file
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	// Create a CSV reader
	reader := csv.NewReader(file)

	// Read the header row
	headers, err := reader.Read()
	if err != nil {
		return nil, err
	}

	// Initialize a list to store maps
	var data []Hit

	// Read the remaining rows
	for {
		row, err := reader.Read()
		if err != nil {
			break // End of file
		}

		// Create a map for the current row
		rowMap := make(Hit)

		// Populate the map with values from the row
		for i, header := range headers {
			if i < len(row) {
				rowMap[header] = row[i]
			}
		}

		// Append the map to the list
		data = append(data, rowMap)
	}

	return data, nil
}

func Blast(query string, name string, db Database, logFile *os.File) ([]Hit, error) {
	var task Task

	switch db.Method {
	case "blastn":
		task = BlastTask{}
	case "diamond":
		task = DiamondTask{}
	case "infernal":
		task = InfernalTask{}
	default:
		return []Hit{}, fmt.Errorf("unknown task: %s", db.Method)
	}

	outFile, err := os.CreateTemp("", "out_*.txt")
	if err != nil {
		return []Hit{}, err
	}
	defer outFile.Close()

	err = task.Run(db.Location, query, outFile.Name(), logFile)
	if err != nil {
		return []Hit{}, err
	}

	err = os.Remove(outFile.Name())
	if err != nil {
		return []Hit{}, err
	}
	return readCSV(outFile.Name())
}

func Annotate(seq string, dbs Databases, linear bool, isDetailed bool) error {
	_, err := getRawHits(seq, linear, dbs)
	return err
}
