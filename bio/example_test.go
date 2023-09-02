package bio_test

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"strings"

	"github.com/TimothyStiles/poly/bio"
	"github.com/TimothyStiles/poly/bio/genbank"
	"github.com/TimothyStiles/poly/bio/gff"
	"github.com/TimothyStiles/poly/bio/polyjson"
)

// This is where the integration tests that make effed up cyclic dependencies go.

func Example() {

	// Poly can take in basic gff, gbk, fasta, and JSON.
	// We call the json package "pson" (poly JSON) to prevent namespace collision with Go's standard json package.

	gffInput, _ := gff.Read("../data/ecoli-mg1655-short.gff")
	gbkInput, _ := genbank.Read("../data/puc19.gbk")
	//fastaInput, _ := fasta.Read("fasta/data/base.fasta")
	jsonInput, _ := polyjson.Read("../data/cat.json")

	// Poly can also output these file formats. Every file format has a corresponding Write function.
	_ = gff.Write(gffInput, "test.gff")
	_ = genbank.Write(gbkInput, "test.gbk")
	//_ = fasta.WriteFile(fastaInput, "test.fasta")
	_ = polyjson.Write(jsonInput, "test.json")

	// Extra tips:

	// 1. All of these file formats can be read and written in JSON format using their native schemas.
	// 2. If you want to convert from one format to another (e.g. genbank to polyjson), you can easily do so with a for-loop and some field mapping.
	// 3. Every file format is unique but they all share a common interface so you can use them with almost every native function in Poly.
}

// ExampleRead shows an example of reading a file from disk.
func ExampleRead() {
	// Read lets you read files from disk into a parser.
	parser, _ := bio.Read(bio.Fasta, "fasta/data/base.fasta")
	records, _ := parser.Parse()

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

func ExampleReadGz() {
	// ReadGz lets you read gzipped files into a parser.
	parser, _ := bio.ReadGz(bio.Fasta, "fasta/data/base.fasta.gz")
	records, _ := parser.Parse()

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

func ExampleReadCompressed() {
	// Sometimes you have files with a different compression algorithms. Here
	// is how you would read a compressed file without ReadGz.

	// We first make a decoderFunc with a compatible function type
	decoderFunc := func(r io.Reader) (io.Reader, error) {
		return gzip.NewReader(r)
	}
	parser, _ := bio.ReadCompressed(bio.Fasta, "fasta/data/base.fasta.gz", decoderFunc)
	records, _ := parser.Parse()

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

func ExampleNewParserGz() {
	// First, lets make a file that is gzip'd, represented by this
	// buffer.
	var file bytes.Buffer
	zw := gzip.NewWriter(&file)
	_, _ = zw.Write([]byte(`>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY

>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA
DIDGDGQVNYEEFVQMMTAK*`))
	zw.Close()

	parser, _ := bio.NewParserGz(bio.Fasta, &file) // Make a parser with the gz file
	records, _ := parser.Parse()                   // Parse all data records from file

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

func ExampleNewParserCompressed() {
	// First, lets make a file that is compressed, represented by this
	// buffer.
	var file bytes.Buffer
	zw := gzip.NewWriter(&file)
	_, _ = zw.Write([]byte(`>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY

>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA
DIDGDGQVNYEEFVQMMTAK*`))
	zw.Close()

	// Next, lets make a decoderFunc with a compatible function type
	decoderFunc := func(r io.Reader) (io.Reader, error) {
		return gzip.NewReader(r)
	}

	parser, _ := bio.NewParserCompressed(bio.Fasta, &file, decoderFunc) // Make a parser with the compressed file
	records, _ := parser.Parse()                                        // Parse all data records from file

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

func ExampleFasta() {
	// The following can be replaced with a any io.Reader. For example,
	// `file, err := os.Open(path)` for file would also work.
	file := strings.NewReader(`>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY

>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA
DIDGDGQVNYEEFVQMMTAK*`)
	parser, _ := bio.NewParser(bio.Fasta, file) // Make a parser with the file
	records, _ := parser.Parse()                // Parse all data records from file

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}
