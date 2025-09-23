#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include "DNASequence.h"
#include "GeneticCode.h"
#include "SequenceAnalyzer.h"
#include "PatternFinder.h"
#include "FastaParser.h"

void showMenu();
void analyzeSequenceFromInput();
void loadFromFile();
void showComplementarySequence(const DNASequence& seq);
void showComposition(const DNASequence& seq);
void translateToProtein(const DNASequence& seq);
void findORFs(const DNASequence& seq);
void completeAnalysis(const DNASequence& seq);
void findPatterns(const DNASequence& seq);
void exportResults(const std::string& results, const std::string& filename);
void runTests();

int main() {
    std::cout << "=== DNA FINDER v1.0 ===" << std::endl;
    std::cout << "Herramienta de Análisis de Secuencias de ADN" << std::endl;
    std::cout << "Desarrollado en C++" << std::endl << std::endl;
    
    int option;
    
    do {
        showMenu();
        std::cout << "> Opción: ";
        std::cin >> option;
        std::cin.ignore(); // Clear input buffer
        
        switch (option) {
            case 1:
                analyzeSequenceFromInput();
                break;
            case 2:
                loadFromFile();
                break;
            case 3:
                runTests();
                break;
            case 0:
                std::cout << "¡Gracias por usar DNA Finder!" << std::endl;
                break;
            default:
                std::cout << "Opción inválida. Intente nuevamente." << std::endl;
        }
        
        if (option != 0) {
            std::cout << "\nPresione Enter para continuar...";
            std::cin.get();
        }
        
    } while (option != 0);
    
    return 0;
}

void showMenu() {
    std::cout << "\n=== MENÚ PRINCIPAL ===" << std::endl;
    std::cout << "1. Analizar secuencia (entrada directa)" << std::endl;
    std::cout << "2. Cargar desde archivo FASTA" << std::endl;
    std::cout << "3. Ejecutar casos de prueba" << std::endl;
    std::cout << "0. Salir" << std::endl;
}

void analyzeSequenceFromInput() {
    std::string sequence;
    std::cout << "\n> Ingrese secuencia de ADN: ";
    std::getline(std::cin, sequence);
    
    if (sequence.empty()) {
        std::cout << "Secuencia vacía." << std::endl;
        return;
    }
    
    DNASequence dna(sequence);
    
    if (!dna.isValid()) {
        std::cout << "Error: La secuencia contiene nucleótidos inválidos." << std::endl;
        std::cout << "Solo se permiten: A, T, C, G y nucleótidos ambiguos (N, R, Y, etc.)" << std::endl;
        return;
    }
    
    int analysisOption;
    
    do {
        std::cout << "\n=== ANÁLISIS DISPONIBLES ===" << std::endl;
        std::cout << "1. Secuencia complementaria" << std::endl;
        std::cout << "2. Análisis de composición (GC%)" << std::endl;
        std::cout << "3. Traducir a proteína" << std::endl;
        std::cout << "4. Buscar ORFs" << std::endl;
        std::cout << "5. Buscar patrones" << std::endl;
        std::cout << "6. Análisis completo" << std::endl;
        std::cout << "7. Exportar resultados" << std::endl;
        std::cout << "0. Volver al menú principal" << std::endl;
        std::cout << "> Opción: ";
        
        std::cin >> analysisOption;
        std::cin.ignore();
        
        switch (analysisOption) {
            case 1:
                showComplementarySequence(dna);
                break;
            case 2:
                showComposition(dna);
                break;
            case 3:
                translateToProtein(dna);
                break;
            case 4:
                findORFs(dna);
                break;
            case 5:
                findPatterns(dna);
                break;
            case 6:
                completeAnalysis(dna);
                break;
            case 7: {
                std::string results = SequenceAnalyzer::generateReport(dna);
                std::string filename;
                std::cout << "Nombre del archivo: ";
                std::getline(std::cin, filename);
                exportResults(results, filename);
                break;
            }
            case 0:
                break;
            default:
                std::cout << "Opción inválida." << std::endl;
        }
        
        if (analysisOption != 0) {
            std::cout << "\nPresione Enter para continuar...";
            std::cin.get();
        }
        
    } while (analysisOption != 0);
}

void loadFromFile() {
    std::string filename;
    std::cout << "\n> Nombre del archivo FASTA: ";
    std::getline(std::cin, filename);
    
    if (!FastaParser::isValidFastaFile(filename)) {
        std::cout << "Error: El archivo no existe o no es un archivo FASTA válido." << std::endl;
        return;
    }
    
    std::vector<FastaSequence> sequences = FastaParser::parseFile(filename);
    
    if (sequences.empty()) {
        std::cout << "Error: No se encontraron secuencias en el archivo." << std::endl;
        return;
    }
    
    std::cout << "Se encontraron " << sequences.size() << " secuencia(s):" << std::endl;
    
    for (size_t i = 0; i < sequences.size(); i++) {
        std::cout << (i+1) << ". " << sequences[i].header 
                  << " (" << sequences[i].sequence.length() << " nt)" << std::endl;
    }
    
    int seqChoice;
    std::cout << "> Seleccione secuencia para analizar: ";
    std::cin >> seqChoice;
    std::cin.ignore();
    
    if (seqChoice < 1 || seqChoice > static_cast<int>(sequences.size())) {
        std::cout << "Selección inválida." << std::endl;
        return;
    }
    
    DNASequence dna(sequences[seqChoice-1].sequence);
    std::cout << "\nAnalizando: " << sequences[seqChoice-1].header << std::endl;
    completeAnalysis(dna);
}

void showComplementarySequence(const DNASequence& seq) {
    std::cout << "\n=== SECUENCIAS COMPLEMENTARIAS ===" << std::endl;
    std::cout << "Original:           " << seq.getSequence() << std::endl;
    std::cout << "Complementaria:     " << seq.getComplement() << std::endl;
    std::cout << "Reversa Compl.:     " << seq.getReverseComplement() << std::endl;
}

void showComposition(const DNASequence& seq) {
    std::cout << "\n=== ANÁLISIS DE COMPOSICIÓN ===" << std::endl;
    
    auto counts = seq.getAllCounts();
    int total = seq.getLength();
    
    for (auto& pair : counts) {
        double percentage = (static_cast<double>(pair.second) / total) * 100;
        std::cout << pair.first << ": " << pair.second 
                  << " (" << std::fixed << std::setprecision(2) << percentage << "%)" << std::endl;
    }
    
    std::cout << "Contenido GC: " << std::fixed << std::setprecision(2) 
              << seq.getGCContent() << "%" << std::endl;
    std::cout << "Peso Molecular: ~" << std::fixed << std::setprecision(0) 
              << seq.getMolecularWeight() << " Da" << std::endl;
}

void translateToProtein(const DNASequence& seq) {
    std::cout << "\n=== TRADUCCIÓN A PROTEÍNA ===" << std::endl;
    
    std::string protein = GeneticCode::translateSequence(seq.getSequence());
    std::cout << "Proteína (frame +1): " << protein << std::endl;
    
    std::cout << "\nDetalle de codones:" << std::endl;
    std::cout << GeneticCode::translateSequenceVerbose(seq.getSequence());
}

void findORFs(const DNASequence& seq) {
    std::cout << "\n=== OPEN READING FRAMES (ORFs) ===" << std::endl;
    
    std::vector<ORF> orfs = SequenceAnalyzer::findORFsAllFrames(seq.getSequence(), 10);
    
    if (orfs.empty()) {
        std::cout << "No se encontraron ORFs de longitud mínima 10 aminoácidos." << std::endl;
        return;
    }
    
    std::cout << "Se encontraron " << orfs.size() << " ORF(s):" << std::endl;
    
    for (size_t i = 0; i < orfs.size() && i < 10; i++) {
        const ORF& orf = orfs[i];
        std::cout << "\nORF " << (i+1) << ":" << std::endl;
        std::cout << "  Frame: " << orf.frame << std::endl;
        std::cout << "  Posición: " << orf.start << "-" << orf.end << std::endl;
        std::cout << "  Longitud: " << orf.length << " aminoácidos" << std::endl;
        std::cout << "  Proteína: " << orf.protein << std::endl;
    }
}

void findPatterns(const DNASequence& seq) {
    std::cout << "\n=== BÚSQUEDA DE PATRONES ===" << std::endl;
    
    int patternOption;
    std::cout << "1. Sitios de restricción comunes" << std::endl;
    std::cout << "2. Buscar patrón personalizado" << std::endl;
    std::cout << "> Opción: ";
    std::cin >> patternOption;
    std::cin.ignore();
    
    if (patternOption == 1) {
        std::vector<PatternMatch> matches = PatternFinder::findRestrictionSites(seq.getSequence());
        
        if (matches.empty()) {
            std::cout << "No se encontraron sitios de restricción." << std::endl;
        } else {
            std::cout << "Sitios de restricción encontrados:" << std::endl;
            for (const auto& match : matches) {
                std::cout << "  " << match.pattern << " en posición " << match.position 
                          << ": " << match.matchedSequence << std::endl;
            }
        }
        
    } else if (patternOption == 2) {
        std::string pattern;
        std::cout << "Ingrese patrón (use N para cualquier nucleótido): ";
        std::getline(std::cin, pattern);
        
        std::vector<PatternMatch> matches = PatternFinder::findPatternWithWildcards(seq.getSequence(), pattern);
        
        if (matches.empty()) {
            std::cout << "No se encontraron coincidencias para el patrón: " << pattern << std::endl;
        } else {
            std::cout << "Coincidencias encontradas:" << std::endl;
            for (const auto& match : matches) {
                std::cout << "  Posición " << match.position 
                          << ": " << match.matchedSequence << std::endl;
            }
        }
    }
}

void completeAnalysis(const DNASequence& seq) {
    std::cout << SequenceAnalyzer::generateReport(seq) << std::endl;
}

void exportResults(const std::string& results, const std::string& filename) {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cout << "Error: No se pudo crear el archivo " << filename << std::endl;
        return;
    }
    
    file << results;
    file.close();
    
    std::cout << "Resultados exportados a: " << filename << std::endl;
}

void runTests() {
    std::cout << "\n=== EJECUTANDO CASOS DE PRUEBA ===" << std::endl;
    
    std::vector<std::string> testSequences = {
        "ATCG",
        "ATGAAATAG", 
        "AAAAAAAAAA",
        "ATCGATCGATCGATCG",
        "GCGCGC"
    };
    
    for (const std::string& testSeq : testSequences) {
        std::cout << "\n--- Probando: " << testSeq << " ---" << std::endl;
        
        DNASequence dna(testSeq);
        
        if (!dna.isValid()) {
            std::cout << "Secuencia inválida" << std::endl;
            continue;
        }
        
        std::cout << "Complementaria: " << dna.getComplement() << std::endl;
        std::cout << "GC%: " << std::fixed << std::setprecision(2) << dna.getGCContent() << "%" << std::endl;
        std::cout << "Traducción: " << GeneticCode::translateSequence(testSeq) << std::endl;
        
        std::vector<ORF> orfs = SequenceAnalyzer::findORFs(testSeq, 1);
        std::cout << "ORFs: " << orfs.size() << std::endl;
    }
    
    std::cout << "\n¡Casos de prueba completados!" << std::endl;
}