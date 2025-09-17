#include "MainWindow.h"
#include <QtCore/QDir>
#include <QtCore/QStandardPaths>
#include <QtCore/QTextStream>
#include <QtCore/QFile>
#include <iostream>
#include <sstream>
#include <iomanip>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , m_centralWidget(nullptr)
    , m_currentSequence(nullptr)
{
    setWindowTitle("DNA Finder v1.0 - Herramienta de Análisis de Secuencias de ADN");
    setMinimumSize(1000, 700);
    resize(1200, 800);
    
    setupMenuBar();
    setupUI();
    setupStatusBar();
    
    updateStatus("Listo. Ingrese una secuencia de ADN para comenzar el análisis.");
}

MainWindow::~MainWindow()
{
    delete m_currentSequence;
}

void MainWindow::setupMenuBar()
{
    // File menu
    QMenu* fileMenu = menuBar()->addMenu("&Archivo");
    
    m_loadAction = new QAction("&Cargar archivo FASTA...", this);
    m_loadAction->setShortcut(QKeySequence::Open);
    connect(m_loadAction, &QAction::triggered, this, &MainWindow::onLoadFromFile);
    fileMenu->addAction(m_loadAction);
    
    m_exportAction = new QAction("&Exportar resultados...", this);
    m_exportAction->setShortcut(QKeySequence::SaveAs);
    connect(m_exportAction, &QAction::triggered, this, &MainWindow::onExportResults);
    fileMenu->addAction(m_exportAction);
    
    fileMenu->addSeparator();
    
    m_exitAction = new QAction("&Salir", this);
    m_exitAction->setShortcut(QKeySequence::Quit);
    connect(m_exitAction, &QAction::triggered, this, &QWidget::close);
    fileMenu->addAction(m_exitAction);
    
    // Help menu
    QMenu* helpMenu = menuBar()->addMenu("&Ayuda");
    
    m_aboutAction = new QAction("&Acerca de DNA Finder...", this);
    connect(m_aboutAction, &QAction::triggered, this, &MainWindow::onAbout);
    helpMenu->addAction(m_aboutAction);
}

void MainWindow::setupStatusBar()
{
    m_statusLabel = new QLabel("Listo");
    m_progressBar = new QProgressBar();
    m_progressBar->setVisible(false);
    
    statusBar()->addWidget(m_statusLabel, 1);
    statusBar()->addWidget(m_progressBar);
}

void MainWindow::setupUI()
{
    m_centralWidget = new QWidget;
    setCentralWidget(m_centralWidget);
    
    QHBoxLayout* mainLayout = new QHBoxLayout(m_centralWidget);
    
    // Create main splitter
    m_mainSplitter = new QSplitter(Qt::Horizontal);
    mainLayout->addWidget(m_mainSplitter);
    
    // Left panel - Input and controls
    QWidget* leftPanel = new QWidget;
    leftPanel->setMaximumWidth(350);
    leftPanel->setMinimumWidth(300);
    
    QVBoxLayout* leftLayout = new QVBoxLayout(leftPanel);
    
    // Input section
    m_inputGroup = new QGroupBox("Entrada de Secuencia");
    QVBoxLayout* inputLayout = new QVBoxLayout(m_inputGroup);
    
    m_sequenceLabel = new QLabel("Ingrese secuencia de ADN:");
    m_sequenceInput = new QTextEdit;
    m_sequenceInput->setMaximumHeight(120);
    m_sequenceInput->setPlaceholderText("Ejemplo: ATGCGATCGTAGCAAGCTG");
    connect(m_sequenceInput, &QTextEdit::textChanged, this, &MainWindow::onSequenceChanged);
    
    m_loadFileButton = new QPushButton("Cargar archivo FASTA");
    connect(m_loadFileButton, &QPushButton::clicked, this, &MainWindow::onLoadFromFile);
    
    inputLayout->addWidget(m_sequenceLabel);
    inputLayout->addWidget(m_sequenceInput);
    inputLayout->addWidget(m_loadFileButton);
    
    leftLayout->addWidget(m_inputGroup);
    
    // Analysis section
    m_analysisGroup = new QGroupBox("Análisis");
    QVBoxLayout* analysisLayout = new QVBoxLayout(m_analysisGroup);
    
    m_analyzeButton = new QPushButton("Análisis Básico");
    m_analyzeButton->setEnabled(false);
    connect(m_analyzeButton, &QPushButton::clicked, this, &MainWindow::onAnalyzeSequence);
    
    m_runAllButton = new QPushButton("Análisis Completo");
    m_runAllButton->setEnabled(false);
    connect(m_runAllButton, &QPushButton::clicked, this, &MainWindow::onRunAllAnalyses);
    
    m_translateButton = new QPushButton("Traducir a Proteína");
    m_translateButton->setEnabled(false);
    connect(m_translateButton, &QPushButton::clicked, this, &MainWindow::onTranslateProtein);
    
    m_findORFsButton = new QPushButton("Buscar ORFs");
    m_findORFsButton->setEnabled(false);
    connect(m_findORFsButton, &QPushButton::clicked, this, &MainWindow::onFindORFs);
    
    m_findPatternsButton = new QPushButton("Buscar Patrones");
    m_findPatternsButton->setEnabled(false);
    connect(m_findPatternsButton, &QPushButton::clicked, this, &MainWindow::onFindPatterns);
    
    m_analyzeCodonsButton = new QPushButton("Análisis de Codones PRO");
    m_analyzeCodonsButton->setEnabled(false);
    m_analyzeCodonsButton->setStyleSheet("QPushButton { background-color: #4CAF50; color: white; font-weight: bold; }");
    connect(m_analyzeCodonsButton, &QPushButton::clicked, this, &MainWindow::onAnalyzeCodons);
    
    analysisLayout->addWidget(m_analyzeButton);
    analysisLayout->addWidget(m_runAllButton);
    analysisLayout->addWidget(m_translateButton);
    analysisLayout->addWidget(m_findORFsButton);
    analysisLayout->addWidget(m_findPatternsButton);
    analysisLayout->addWidget(m_analyzeCodonsButton);
    analysisLayout->addStretch();
    
    leftLayout->addWidget(m_analysisGroup);
    leftLayout->addStretch();
    
    m_mainSplitter->addWidget(leftPanel);
    
    // Right panel - Results
    QWidget* rightPanel = new QWidget;
    QVBoxLayout* rightLayout = new QVBoxLayout(rightPanel);
    
    QLabel* resultsLabel = new QLabel("Resultados del Análisis:");
    resultsLabel->setStyleSheet("font-weight: bold; font-size: 14px;");
    
    m_resultsTab = new QTabWidget;
    
    // Basic info tab
    m_basicInfoText = new QTextEdit;
    m_basicInfoText->setReadOnly(true);
    m_resultsTab->addTab(m_basicInfoText, "Información Básica");
    
    // Composition tab
    m_compositionText = new QTextEdit;
    m_compositionText->setReadOnly(true);
    m_resultsTab->addTab(m_compositionText, "Composición");
    
    // Translation tab
    m_translationText = new QTextEdit;
    m_translationText->setReadOnly(true);
    m_resultsTab->addTab(m_translationText, "Traducción");
    
    // ORFs tab
    m_orfsText = new QTextEdit;
    m_orfsText->setReadOnly(true);
    m_resultsTab->addTab(m_orfsText, "ORFs");
    
    // Patterns tab
    m_patternsText = new QTextEdit;
    m_patternsText->setReadOnly(true);
    m_resultsTab->addTab(m_patternsText, "Patrones");
    
    // Codon Analysis tab
    m_codonAnalysisText = new QTextEdit;
    m_codonAnalysisText->setReadOnly(true);
    m_codonAnalysisText->setStyleSheet("QTextEdit { font-family: 'Courier New', monospace; }");
    m_resultsTab->addTab(m_codonAnalysisText, "🔬 Análisis Codones PRO");
    
    // Complete report tab
    m_completeReportText = new QTextEdit;
    m_completeReportText->setReadOnly(true);
    m_resultsTab->addTab(m_completeReportText, "Reporte Completo");
    
    rightLayout->addWidget(resultsLabel);
    rightLayout->addWidget(m_resultsTab);
    
    m_mainSplitter->addWidget(rightPanel);
    m_mainSplitter->setSizes({350, 650});
}

void MainWindow::onSequenceChanged()
{
    QString text = m_sequenceInput->toPlainText().trimmed().toUpper();
    bool hasSequence = !text.isEmpty();
    
    m_analyzeButton->setEnabled(hasSequence);
    m_runAllButton->setEnabled(hasSequence);
    m_translateButton->setEnabled(hasSequence);
    m_findORFsButton->setEnabled(hasSequence);
    m_findPatternsButton->setEnabled(hasSequence);
    m_analyzeCodonsButton->setEnabled(hasSequence);
    
    if (hasSequence) {
        updateStatus(QString("Secuencia ingresada: %1 nucleótidos").arg(text.length()));
    } else {
        updateStatus("Ingrese una secuencia de ADN para comenzar el análisis.");
        clearResults();
    }
}

void MainWindow::onAnalyzeSequence()
{
    QString sequence = m_sequenceInput->toPlainText().trimmed().toUpper();
    
    if (sequence.isEmpty()) {
        QMessageBox::warning(this, "Advertencia", "Por favor ingrese una secuencia de ADN.");
        return;
    }
    
    delete m_currentSequence;
    m_currentSequence = new DNASequence(sequence.toStdString());
    
    if (!m_currentSequence->isValid()) {
        QMessageBox::critical(this, "Error", 
            "La secuencia contiene nucleótidos inválidos.\n"
            "Solo se permiten: A, T, C, G y nucleótidos ambiguos (N, R, Y, etc.)");
        return;
    }
    
    displaySequenceInfo(*m_currentSequence);
    displayComposition(*m_currentSequence);
    
    m_resultsTab->setCurrentIndex(0);
    updateStatus("Análisis básico completado.");
}

void MainWindow::onRunAllAnalyses()
{
    QString sequence = m_sequenceInput->toPlainText().trimmed().toUpper();
    
    if (sequence.isEmpty()) {
        QMessageBox::warning(this, "Advertencia", "Por favor ingrese una secuencia de ADN.");
        return;
    }
    
    m_progressBar->setVisible(true);
    m_progressBar->setValue(0);
    
    delete m_currentSequence;
    m_currentSequence = new DNASequence(sequence.toStdString());
    
    if (!m_currentSequence->isValid()) {
        QMessageBox::critical(this, "Error", 
            "La secuencia contiene nucleótidos inválidos.\n"
            "Solo se permiten: A, T, C, G y nucleótidos ambiguos (N, R, Y, etc.)");
        m_progressBar->setVisible(false);
        return;
    }
    
    updateStatus("Ejecutando análisis completo...");
    
    displaySequenceInfo(*m_currentSequence);
    m_progressBar->setValue(20);
    
    displayComposition(*m_currentSequence);
    m_progressBar->setValue(40);
    
    displayTranslation(*m_currentSequence);
    m_progressBar->setValue(60);
    
    displayORFs(*m_currentSequence);
    m_progressBar->setValue(70);
    
    displayPatterns(*m_currentSequence);
    m_progressBar->setValue(75);
    
    displayCodonAnalysis(*m_currentSequence);
    m_progressBar->setValue(90);
    
    // Generate complete report
    std::string report = SequenceAnalyzer::generateReport(*m_currentSequence);
    m_completeReportText->setPlainText(QString::fromStdString(report));
    m_currentResults = QString::fromStdString(report);
    
    m_progressBar->setValue(100);
    m_progressBar->setVisible(false);
    
    m_resultsTab->setCurrentIndex(6); // Show complete report
    updateStatus("Análisis completo finalizado.");
}

void MainWindow::onTranslateProtein()
{
    if (!m_currentSequence) {
        onAnalyzeSequence();
        if (!m_currentSequence) return;
    }
    
    displayTranslation(*m_currentSequence);
    m_resultsTab->setCurrentIndex(2);
    updateStatus("Traducción a proteína completada.");
}

void MainWindow::onFindORFs()
{
    if (!m_currentSequence) {
        onAnalyzeSequence();
        if (!m_currentSequence) return;
    }
    
    displayORFs(*m_currentSequence);
    m_resultsTab->setCurrentIndex(3);
    updateStatus("Búsqueda de ORFs completada.");
}

void MainWindow::onFindPatterns()
{
    if (!m_currentSequence) {
        onAnalyzeSequence();
        if (!m_currentSequence) return;
    }
    
    displayPatterns(*m_currentSequence);
    m_resultsTab->setCurrentIndex(4);
    updateStatus("Búsqueda de patrones completada.");
}

void MainWindow::onAnalyzeCodons()
{
    if (!m_currentSequence) {
        onAnalyzeSequence();
        if (!m_currentSequence) return;
    }
    
    displayCodonAnalysis(*m_currentSequence);
    m_resultsTab->setCurrentIndex(5);  // Codon analysis tab
    updateStatus("Análisis profesional de codones completado.");
}

void MainWindow::onLoadFromFile()
{
    QString fileName = QFileDialog::getOpenFileName(this,
        "Cargar archivo FASTA", 
        QDir::currentPath(),
        "Archivos FASTA (*.fasta *.fas *.fa);;Todos los archivos (*.*)");
    
    if (fileName.isEmpty()) {
        return;
    }
    
    if (!FastaParser::isValidFastaFile(fileName.toStdString())) {
        QMessageBox::critical(this, "Error", 
            "El archivo no existe o no es un archivo FASTA válido.");
        return;
    }
    
    std::vector<FastaSequence> sequences = FastaParser::parseFile(fileName.toStdString());
    
    if (sequences.empty()) {
        QMessageBox::warning(this, "Advertencia", 
            "No se encontraron secuencias en el archivo.");
        return;
    }
    
    if (sequences.size() == 1) {
        m_sequenceInput->setPlainText(QString::fromStdString(sequences[0].sequence));
        updateStatus(QString("Cargada secuencia: %1").arg(QString::fromStdString(sequences[0].header)));
    } else {
        // Multiple sequences - show first one and inform user
        m_sequenceInput->setPlainText(QString::fromStdString(sequences[0].sequence));
        QMessageBox::information(this, "Información", 
            QString("Se encontraron %1 secuencias. Se cargó la primera: %2")
            .arg(sequences.size())
            .arg(QString::fromStdString(sequences[0].header)));
        updateStatus(QString("Cargadas %1 secuencias del archivo").arg(sequences.size()));
    }
}

void MainWindow::onExportResults()
{
    if (m_currentResults.isEmpty()) {
        QMessageBox::warning(this, "Advertencia", 
            "No hay resultados para exportar. Ejecute un análisis primero.");
        return;
    }
    
    QString fileName = QFileDialog::getSaveFileName(this,
        "Exportar resultados", 
        "dna_analysis_results.txt",
        "Archivos de texto (*.txt);;Todos los archivos (*.*)");
    
    if (fileName.isEmpty()) {
        return;
    }
    
    QFile file(fileName);
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream stream(&file);
        stream << m_currentResults;
        file.close();
        
        QMessageBox::information(this, "Éxito", 
            QString("Resultados exportados exitosamente a:\n%1").arg(fileName));
        updateStatus(QString("Resultados exportados a: %1").arg(fileName));
    } else {
        QMessageBox::critical(this, "Error", 
            QString("No se pudo escribir el archivo:\n%1").arg(fileName));
    }
}

void MainWindow::onAbout()
{
    QMessageBox::about(this, "Acerca de DNA Finder",
        "<h2>DNA Finder v1.0</h2>"
        "<p>Herramienta de Análisis de Secuencias de ADN</p>"
        "<p>Desarrollado en C++ con Qt</p>"
        "<p><b>Funcionalidades:</b></p>"
        "<ul>"
        "<li>Validación de secuencias de ADN</li>"
        "<li>Cálculo de secuencias complementarias</li>"
        "<li>Análisis de composición nucleotídica</li>"
        "<li>Traducción a secuencias de aminoácidos</li>"
        "<li>Detección de marcos de lectura abiertos (ORFs)</li>"
        "<li>Búsqueda de patrones y sitios de restricción</li>"
        "</ul>"
        "<p>Ideal para estudiantes y profesionales en biotecnología.</p>");
}

void MainWindow::displaySequenceInfo(const DNASequence& seq)
{
    QString info;
    info += QString("=== INFORMACIÓN BÁSICA ===\n\n");
    info += QString("Secuencia original: %1\n").arg(QString::fromStdString(seq.getSequence()));
    info += QString("Longitud: %1 nucleótidos\n\n").arg(seq.getLength());
    info += QString("Secuencia complementaria: %1\n").arg(QString::fromStdString(seq.getComplement()));
    info += QString("Reversa complementaria: %1\n\n").arg(QString::fromStdString(seq.getReverseComplement()));
    info += QString("Validez: %1\n").arg(seq.isValid() ? "✓ Válida" : "✗ Inválida");
    
    m_basicInfoText->setPlainText(info);
}

void MainWindow::displayComposition(const DNASequence& seq)
{
    QString composition;
    composition += QString("=== ANÁLISIS DE COMPOSICIÓN ===\n\n");
    
    auto counts = seq.getAllCounts();
    int total = seq.getLength();
    
    composition += QString("Frecuencia de nucleótidos:\n");
    for (auto& pair : counts) {
        double percentage = (static_cast<double>(pair.second) / total) * 100;
        composition += QString("• %1: %2 (%3%)\n")
                      .arg(pair.first)
                      .arg(pair.second)
                      .arg(percentage, 0, 'f', 2);
    }
    
    composition += QString("\nContenido GC: %1%\n").arg(seq.getGCContent(), 0, 'f', 2);
    composition += QString("Peso Molecular: ~%1 Da\n").arg(seq.getMolecularWeight(), 0, 'f', 0);
    
    m_compositionText->setPlainText(composition);
}

void MainWindow::displayTranslation(const DNASequence& seq)
{
    QString translation;
    translation += QString("=== TRADUCCIÓN A PROTEÍNA ===\n\n");
    
    std::string protein = GeneticCode::translateSequence(seq.getSequence());
    translation += QString("Secuencia de aminoácidos: %1\n\n").arg(QString::fromStdString(protein));
    
    translation += QString("Detalle de codones:\n");
    translation += QString::fromStdString(GeneticCode::translateSequenceVerbose(seq.getSequence()));
    
    m_translationText->setPlainText(translation);
}

void MainWindow::displayORFs(const DNASequence& seq)
{
    QString orfs;
    orfs += QString("=== MARCOS DE LECTURA ABIERTOS (ORFs) ===\n\n");
    
    std::vector<ORF> orfList = SequenceAnalyzer::findORFsAllFrames(seq.getSequence(), 10);
    
    if (orfList.empty()) {
        orfs += QString("No se encontraron ORFs de longitud mínima 10 aminoácidos.\n");
    } else {
        orfs += QString("Se encontraron %1 ORF(s):\n\n").arg(orfList.size());
        
        for (size_t i = 0; i < orfList.size() && i < 20; i++) {
            const ORF& orf = orfList[i];
            orfs += QString("ORF %1:\n").arg(i + 1);
            orfs += QString("  • Frame: %1\n").arg(orf.frame);
            orfs += QString("  • Posición: %1-%2\n").arg(orf.start).arg(orf.end);
            orfs += QString("  • Longitud: %1 aminoácidos\n").arg(orf.length);
            orfs += QString("  • Proteína: %1\n\n").arg(QString::fromStdString(orf.protein));
        }
        
        if (orfList.size() > 20) {
            orfs += QString("... y %1 ORFs adicionales\n").arg(orfList.size() - 20);
        }
    }
    
    m_orfsText->setPlainText(orfs);
}

void MainWindow::displayPatterns(const DNASequence& seq)
{
    QString patterns;
    patterns += QString("=== BÚSQUEDA DE PATRONES ===\n\n");
    
    // Search for restriction sites
    std::vector<PatternMatch> matches = PatternFinder::findRestrictionSites(seq.getSequence());
    
    patterns += QString("Sitios de restricción encontrados:\n");
    if (matches.empty()) {
        patterns += QString("• No se encontraron sitios de restricción comunes.\n\n");
    } else {
        for (const auto& match : matches) {
            patterns += QString("• %1 en posición %2: %3\n")
                       .arg(QString::fromStdString(match.pattern))
                       .arg(match.position)
                       .arg(QString::fromStdString(match.matchedSequence));
        }
        patterns += QString("\n");
    }
    
    // Additional pattern analysis could be added here
    patterns += QString("Para búsquedas de patrones personalizados,\n");
    patterns += QString("use la versión de línea de comandos del programa.\n");
    
    m_patternsText->setPlainText(patterns);
}

void MainWindow::displayCodonAnalysis(const DNASequence& seq)
{
    CodonAnalyzer analyzer;
    CodonAnalysisReport report = analyzer.analyzeCodonUsage(seq.getSequence(), "E.coli");
    
    QString analysis = QString::fromStdString(analyzer.generateCodonReport(report));
    m_codonAnalysisText->setPlainText(analysis);
}

void MainWindow::clearResults()
{
    m_basicInfoText->clear();
    m_compositionText->clear();
    m_translationText->clear();
    m_orfsText->clear();
    m_patternsText->clear();
    m_codonAnalysisText->clear();
    m_completeReportText->clear();
    m_currentResults.clear();
}

void MainWindow::updateStatus(const QString& message)
{
    m_statusLabel->setText(message);
}