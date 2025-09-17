#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QLabel>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QSplitter>
#include <QtCore/QTimer>

#include "DNASequence.h"
#include "SequenceAnalyzer.h"
#include "GeneticCode.h"
#include "PatternFinder.h"
#include "FastaParser.h"
#include "CodonAnalyzer.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void onAnalyzeSequence();
    void onLoadFromFile();
    void onExportResults();
    void onAbout();
    void onSequenceChanged();
    void onRunAllAnalyses();
    void onFindORFs();
    void onFindPatterns();
    void onTranslateProtein();
    void onAnalyzeCodons();

private:
    void setupUI();
    void setupMenuBar();
    void setupStatusBar();
    void setupCentralWidget();
    
    void displaySequenceInfo(const DNASequence& seq);
    void displayComposition(const DNASequence& seq);
    void displayTranslation(const DNASequence& seq);
    void displayORFs(const DNASequence& seq);
    void displayPatterns(const DNASequence& seq);
    void displayCodonAnalysis(const DNASequence& seq);
    
    void clearResults();
    void updateStatus(const QString& message);
    
    // UI Components
    QWidget* m_centralWidget;
    QSplitter* m_mainSplitter;
    
    // Input section
    QGroupBox* m_inputGroup;
    QTextEdit* m_sequenceInput;
    QPushButton* m_loadFileButton;
    QPushButton* m_analyzeButton;
    QLabel* m_sequenceLabel;
    
    // Analysis buttons
    QGroupBox* m_analysisGroup;
    QPushButton* m_runAllButton;
    QPushButton* m_findORFsButton;
    QPushButton* m_translateButton;
    QPushButton* m_findPatternsButton;
    QPushButton* m_analyzeCodonsButton;
    
    // Results section
    QTabWidget* m_resultsTab;
    QTextEdit* m_basicInfoText;
    QTextEdit* m_compositionText;
    QTextEdit* m_translationText;
    QTextEdit* m_orfsText;
    QTextEdit* m_patternsText;
    QTextEdit* m_codonAnalysisText;
    QTextEdit* m_completeReportText;
    
    // Status and progress
    QProgressBar* m_progressBar;
    QLabel* m_statusLabel;
    
    // Menu and actions
    QAction* m_loadAction;
    QAction* m_exportAction;
    QAction* m_exitAction;
    QAction* m_aboutAction;
    
    // Current sequence
    DNASequence* m_currentSequence;
    QString m_currentResults;
};

#endif // MAINWINDOW_H