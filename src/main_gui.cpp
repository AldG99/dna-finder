#include <QtWidgets/QApplication>
#include <QtCore/QDir>
#include <QtCore/QStandardPaths>
#include <iostream>

#include "MainWindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    
    // Set application properties
    app.setApplicationName("DNA Finder");
    app.setApplicationVersion("1.0");
    app.setOrganizationName("DNA Finder Development");
    app.setApplicationDisplayName("DNA Finder - Herramienta de Análisis de Secuencias de ADN");
    
    // Set the working directory to where the executable is located
    QDir::setCurrent(QCoreApplication::applicationDirPath());
    
    try {
        MainWindow window;
        window.show();
        
        return app.exec();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
}