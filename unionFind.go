package main

import (
	"bufio"
	"errors"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

const centerAtom = 6285

func getDistances() []float64 {
	return []float64{10, 20, 30, 40, 50, 60}
}

func createStructuresDir(dir string) {

	fileInfo, err := ioutil.ReadDir(dir)
	if err != nil {
		fmt.Println("failed to read directory: " + dir)
		log.Fatal(err)
	}
	for _, file := range fileInfo {
		if !strings.Contains(file.Name(), "dist") && !strings.Contains(file.Name(), "ion") {
			fmt.Println("Creating derivative structures from " + file.Name())
			createStructures(dir, file.Name())
		}
	}
}

func createStructures(dir string, file string) {
	dists := getDistances()
	createStructure(dir, file, -1, true)
	for _, dist := range dists {
		createStructure(dir, file, dist, true)
		createStructure(dir, file, dist, false)
	}
}

func createStructure(dir string, inName string, distance float64, deleteIon bool) {

	frameName, specs, atoms := loadFrame(filepath.Join(dir, inName))
	introPrint := "Creating structure from " + inName
	if deleteIon {
		introPrint += " with no ion"
	}
	if distance > 0 {
		introPrint += " and with a cutoff dist of " + strconv.Itoa(int(distance))
	}
	fmt.Println(introPrint)
	if deleteIon {
		atoms[centerAtom].markedForDelete = true
	}
	//centerPos := atoms[centerAtom].pos

	if distance > 0.0 {
		fmt.Println("Marking atoms outside new box size")
		markOutsideAtoms(atoms, distance)
	}
	fmt.Println("Identifying molecules via union-find")
	consolidateMolecules(atoms)
	fmt.Println("Deleting atoms that are part of molecules with any atom outside new box")
	deleteAtomsNotOfGroup(atoms)
	/*if distance > 0.0 {
		fmt.Println("Shifting center")
		shiftCenter(atoms, centerPos)
		//specs = getSpecs(atoms)
	}*/
	fmt.Println("Renumbering atoms")
	newAtoms := copyMoleculeWithRenumbering(atoms)

	outName := frameName
	if distance > 0 {
		outName += "_dist" + strconv.Itoa(int(distance))
	}
	if deleteIon {
		outName += "_noion"
	}
	writeFragment(newAtoms, dir, outName, specs)
}

/*func getSpecs(atoms []*atom) float64 {
	minX := 0
	maxX := 0
	minY := 0
	maxY := 0
	minZ := 0
	maxZ := 0

}*/

func writeFragment(atoms []*atom, dir string, name string, specs string) {
	os.MkdirAll(dir, 0755)
	thisPath := filepath.Join(dir, name + ".txyz")
	// fmt.Println(thisPath)
	thisFile, err := os.Create(thisPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + thisPath)
		log.Fatal(err)
	}

	// write header
	_, err = thisFile.WriteString(strconv.Itoa(len(atoms)) + "\t Fragment " + name + "\n")
	if err != nil {
		fmt.Println("Failed to write header line to key: " + thisPath)
		log.Fatal(err)
	}
	_, err = thisFile.WriteString(specs + "\n")
	// write body
	for i := 0; i < len(atoms); i++ {
		line := strconv.Itoa(i+1) + "\t" + atoms[i].element + "\t" + fmt.Sprintf("%.6f",atoms[i].pos[0]) + "\t" +
			fmt.Sprintf("%.6f",atoms[i].pos[1]) + "\t" + fmt.Sprintf("%.6f",atoms[i].pos[2]) + "\t" +
			strconv.Itoa(atoms[i].atomType)
		for _, bondedAtom := range atoms[i].bondedAtoms {
			line += "\t" + strconv.Itoa(bondedAtom)
		}

		_, err = thisFile.WriteString(line + "\n")
		if err != nil {
			fmt.Println("Failed to write header line to key: " + thisPath)
			log.Fatal(err)
		}
	}

}

func shiftCenter(atoms map[int]*atom, centerPos []float64) {
	for _, atom := range atoms {
		atom.pos[0] = atom.pos[0] - centerPos[0]
		atom.pos[1] = atom.pos[1] - centerPos[1]
		atom.pos[2] = atom.pos[2] - centerPos[2]
	}
}

func consolidateMolecules(atoms map[int]*atom) {
	for atomID, atom := range atoms {
		for _, bondedAtom := range atom.bondedAtoms {
			union(atoms, atomID, bondedAtom)
		}
	}

}

func markOutsideAtoms(atoms map[int]*atom, distance float64) {
	center := atoms[6285].pos
	for _, atom := range atoms {
		/*if math.Abs(atom.pos[0]-center[0]) > distance || math.Abs(atom.pos[1]-center[1]) > distance || math.Abs(atom.pos[2]-center[2]) > distance {
			atom.markedForDelete = true
		}*/
		if getDistance(center, atom.pos) > distance && isAtomDeletionEligible(atom.atomType) {
			atom.markedForDelete = true
		}
	}
}

func isAtomDeletionEligible(i int) bool {
	if i < 800 && i >= 700 { // lipid
		return true
	} else if i == 349 || i == 350 { // water
		return true
	} else if i == 361 || i == 353 { // Cl- or K+ ions
		return true
	} else {
		return false
	}
}

func getDistance(p1 []float64, p2 []float64) float64 {
	return math.Sqrt(math.Pow(p1[0]-p2[0],2) + math.Pow(p1[1]-p2[1],2) + math.Pow(p1[2]-p2[2],2))
}

func loadFrame(filePath string) (string, string, map[int]*atom) {

	// Create structure to store atoms
	atoms := make(map[int]*atom)
	var specs string
	// open file
	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("Failed to open molecule file: " + filePath)
		log.Fatal(err)
	}
	lipidName := strings.Split(filepath.Base(filePath),".")[0]


	// Initialize scanner
	scanner := bufio.NewScanner(file)
	// ignore first line
	scanner.Scan()
	// create line counter
	i := 1
	// iterate over all other lines
	for scanner.Scan() {
		// get next line
		line := scanner.Text()
		// split by whitespace
		tokens := strings.Fields(line)
		// check line length before proceeding
		if len(tokens) >= 6 && !strings.Contains(line, "90.000000   90.000000   90.000000") {

			// create new atom
			var newAtom atom

			// get number of atom from file
			atomNum, err := strconv.Atoi(tokens[0])
			if err != nil {
				newErr := errors.New("Failed to convert token in position 0 on line " + strconv.Itoa(i) + " to an integer")
				log.Fatal(newErr)
			}

			// assign parent and size for union-find algorithm
			newAtom.id = atomNum
			newAtom.parent = atomNum
			newAtom.treeSize = 1

			// assign element
			newAtom.element = tokens[1]

			// assign positions
			pos := make([]float64,3)
			for j := 2; j < 5; j++ {
				pos[j-2], err = strconv.ParseFloat(tokens[j],64)
				if err != nil {
					newErr := errors.New("Failed to convert token in position 0 on line " + strconv.Itoa(j) + " to a float64")
					log.Fatal(newErr)
				}
			}
			newAtom.pos = pos

			// assign atomType from file
			newAtom.atomType, err = strconv.Atoi(tokens[5])
			if err != nil {
				newErr := errors.New("Failed to convert token in position 5 on line " + strconv.Itoa(i) + " to an integer")
				log.Fatal(newErr)
			}

			// assign bonds from file
			bonds := make([]int,len(tokens)-6)
			for j := 6; j < len(tokens); j++ {
				bonds[j-6], err = strconv.Atoi(tokens[j])
				if err != nil {
					newErr := errors.New("Failed to convert token in position " + strconv.Itoa(j) + " on line " + strconv.Itoa(i) + " to an integer")
					log.Fatal(newErr)
				}
			}
			newAtom.bondedAtoms = bonds

			// add atom to map
			atoms[atomNum] = &newAtom

		} else if len(tokens) >= 6 && strings.Contains(line, "90.000000   90.000000   90.000000") {
			specs = line
		}
		i++
	}
	fmt.Println("Read in lipid with " + strconv.Itoa(len(atoms)) + " atoms")
	return lipidName, specs, atoms
}

////////////////
// Union Find Alg
////////////////

type atom struct {
	element  string
	pos []float64
	atomType int
	id       int
	bondedAtoms []int
	parent int
	treeSize int
	markedForDelete bool
}

func union(atoms map[int]*atom, atom1 int, atom2 int) {
	root1 := root(atoms, atom1)
	root2 := root(atoms, atom2)
	if root1 != root2 {
		if atoms[root1].treeSize < atoms[root2].treeSize {
			atoms[root1].parent = root2
			atoms[root2].treeSize += atoms[root1].treeSize
		} else {
			atoms[root2].parent = root1
			atoms[root1].treeSize += atoms[atom1].treeSize
		}
	}
}

func connected(atoms map[int]*atom, atom1 int, atom2 int) bool {
	root1 := root(atoms, atom1)
	root2 := root(atoms, atom2)
	if root1 != root2 {
		return false
	} else {
		return true
	}
}

func validate(atoms map[int]*atom, atom1 int) bool {
	if _, ok := atoms[atom1]; ok {
		return true
	}
	return false
}

func root(atoms map[int]*atom, atom1 int) int {
	isAtomInMap := validate(atoms, atom1)

	if isAtomInMap {
		// store array of visited atoms for path compression afterwards
		var visitedAtoms []int

		// check if atom's parent is equal to atom's parent's parent (i.e. we have reached the top of the tree)
		for atoms[atom1].parent != atoms[atoms[atom1].parent].parent {
			// if not set atom's parent to atom's parent's parent
			atoms[atom1].parent = atoms[atoms[atom1].parent].parent
			// save atom's parent to list of visited atoms to path compress afterwards
			visitedAtoms = append(visitedAtoms, atoms[atom1].parent)
		}

		// compress path
		for _, visitedAtom := range visitedAtoms {
			atoms[visitedAtom].parent = atoms[atom1].parent
		}

		return atoms[atom1].parent
	} else {
		err := errors.New("Call to root(): atom ID key " + strconv.Itoa(atom1) + " was not found in molecule map.")
		log.Fatal(err)
		return -1
	}
}

// Includes renumbering atom map keys from 1 to len(molecule)-1
func copyMoleculeWithRenumbering(atomsA map[int]*atom) []*atom {
	var atomsSlice []*atom
	for _, atom := range atomsA {
		atomsSlice = append(atomsSlice, atom)
	}
	qsort(atomsSlice)
	atomIDOldToNewMap := make(map[string]string)

	for i, atom := range atomsSlice {
		atomIDOldToNewMap[strconv.Itoa(atom.id)] = strconv.Itoa(i + 1)
		i++
	}

	for _, thisAtom := range atomsSlice {
		for j := 0; j < len(thisAtom.bondedAtoms); j++ {
			newBondID, _ := strconv.Atoi(atomIDOldToNewMap[strconv.Itoa(thisAtom.bondedAtoms[j])])
			thisAtom.bondedAtoms[j] = newBondID
		}
	}
	return atomsSlice
}

func deleteAtomsNotOfGroup(atoms map[int]*atom) {
	for atomID, atom := range atoms {
		if atom.markedForDelete {
			root1 := root(atoms, atomID)
			atoms[root1].markedForDelete = true
		}
	}
	// delete all atoms not of the root of given groups
	atomIDs := make([]int, len(atoms))
	markForDelete := make([]bool, len(atoms))

	i := 0
	deletions := 0
	for atomID := range atoms {
		atomIDs[i] = atomID
		if atoms[root(atoms, atomID)].markedForDelete {
			markForDelete[i] = true
			deletions++
		}
		i++
	}
	for k := 0; k < len(atomIDs); k++ {
		if markForDelete[k] == true {
			delete(atoms, atomIDs[k])
		}
	}
	fmt.Println("Deleted " + strconv.Itoa(deletions) + " atoms")
}