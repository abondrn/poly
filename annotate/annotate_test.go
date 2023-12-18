package annotate_test

import (
	"fmt"
	"os"
	"testing"

	"github.com/TimothyStiles/poly/annotate"
	"github.com/stretchr/testify/assert"
)

func TestLoadResources(t *testing.T) {
	_, err := annotate.LoadDatabases("../../pLannotate/plannotate/data/data/databases.yml")
	assert.NoError(t, err)
}

func RRNB() (string, error) {
	content, err := os.ReadFile("../../pLannotate/tests/test_data/RRNB_fragment.txt")
	if err != nil {
		return "", err
	}
	return string(content), nil
}

func TestBlast(t *testing.T) {
	dbs, _ := annotate.LoadDatabases("../../pLannotate/plannotate/data/data/databases.yml")
	fmt.Printf("%v\n", dbs)

	rrnb, err := RRNB()
	assert.NoError(t, err)

	snapgene_db, exists := dbs["snapgene"]
	assert.True(t, exists, "snapgene database config not found")

	inPath, err := annotate.CreateTempFasta(rrnb)
	assert.NoError(t, err)

	logFile, err := os.Create("output.log")
	assert.NoError(t, err)
	defer logFile.Close()

	hits, err := annotate.Blast(inPath, "snapgene", snapgene_db, logFile)
	logs, err := os.ReadFile("output.log")
	fmt.Println(string(logs))
	assert.NoError(t, err)
	assert.NotEmpty(t, hits)
}

func TestAnnotate(t *testing.T) {
	dbs, _ := annotate.LoadDatabases("../../pLannotate/plannotate/data/data/databases.yml")

	rrnb, err := RRNB()

	err = annotate.Annotate(rrnb, dbs, false, false)
	logs, err := os.ReadFile("output.log")
	fmt.Println(string(logs))
	assert.NoError(t, err)
	//assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"

	err = annotate.Annotate(rrnb, dbs, true, false)
	logs, err = os.ReadFile("output.log")
	fmt.Println(string(logs))
	assert.NoError(t, err)
	//assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"
}
