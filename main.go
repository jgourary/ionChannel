package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

func main() {
	traj2frames("input\\kcsa.arc", "output\\frames")
	// createStructuresDir("output\\frames")
	// createStructures("", "frame_10.txyz")
}

func traj2frames(filePath string, outPath string) {
	// open file
	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("Failed to open molecule file: " + filePath)
		log.Fatal(err)
	}

	scanner := bufio.NewScanner(file)

	var atoms []string
	specs := ""
	frameCount := 0
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Fields(line)
		if len(fields) > 5 && !strings.Contains(line, "90.000000   90.000000   90.000000") {
			atoms = append(atoms, line)
		} else if len(fields) > 5 && strings.Contains(line, "90.000000   90.000000   90.000000") {
			specs = line
		} else if len(fields) == 1 {
			frameCount++
			fmt.Println("Processed frame " + strconv.Itoa(frameCount))
			atoms2txyz(atoms, specs, outPath, frameCount)
			atoms = []string{}

		}
	}
}

func atoms2txyz(atoms []string, specs string, outPath string, frameCount int) {
	outFile, err := os.Create(filepath.Join(outPath, "frame_"+strconv.Itoa(frameCount)+".txyz"))
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + outPath)
		log.Fatal(err)
	}
	_, _ = outFile.WriteString(strconv.Itoa(len(atoms)) + "\n")
	_, _ = outFile.WriteString(specs + "\n")
	for _, atom := range atoms {
		_, _ = outFile.WriteString(atom + "\n")
	}
}
