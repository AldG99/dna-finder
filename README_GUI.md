# DNA Finder - Versión con Interfaz Gráfica

![DNA Finder GUI](https://img.shields.io/badge/GUI-Qt5/Qt6-green) ![C++](https://img.shields.io/badge/C%2B%2B-11%2B-blue) ![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey)

## 🖥️ Descripción

Esta es la versión con interfaz gráfica del DNA Finder, que proporciona una experiencia visual intuitiva para el análisis de secuencias de ADN.

## ✨ Características de la GUI

### Interfaz Principal
- **Panel de entrada**: Área de texto para escribir secuencias de ADN
- **Carga de archivos**: Botón para cargar archivos FASTA
- **Panel de análisis**: Botones para diferentes tipos de análisis
- **Resultados por pestañas**: Organización clara de los resultados

### Funcionalidades
- ✅ Entrada directa de secuencias
- ✅ Carga de archivos FASTA
- ✅ Validación en tiempo real
- ✅ Análisis básico y completo
- ✅ Traducción a proteínas
- ✅ Búsqueda de ORFs
- ✅ Detección de patrones
- ✅ Exportación de resultados
- ✅ Barra de progreso para análisis largos

## 📋 Requisitos del Sistema

### Dependencias Necesarias

#### macOS
```bash
# Instalar Qt usando Homebrew
brew install qt

# O instalar Qt6 específicamente
brew install qt6
```

#### Ubuntu/Debian
```bash
# Qt5 (recomendado)
sudo apt update
sudo apt install qtbase5-dev qtbase5-dev-tools

# O Qt6
sudo apt install qt6-base-dev qt6-base-dev-tools
```

#### Fedora/RHEL
```bash
# Qt5
sudo dnf install qt5-qtbase-devel

# O Qt6  
sudo dnf install qt6-qtbase-devel
```

#### Arch Linux
```bash
# Qt5
sudo pacman -S qt5-base

# O Qt6
sudo pacman -S qt6-base
```

#### Windows
1. Descargar Qt desde [qt.io](https://www.qt.io/download)
2. Instalar Qt Creator con las librerías de desarrollo
3. Configurar el PATH para incluir las herramientas de Qt

### Herramientas Adicionales
```bash
# Verificar que pkg-config esté instalado
pkg-config --version

# En Ubuntu/Debian si no está instalado
sudo apt install pkg-config
```

## 🔧 Compilación

### 1. Verificar Dependencias
```bash
make -f Makefile.gui check-qt
```

### 2. Compilar la Versión GUI
```bash
# Solo GUI
make -f Makefile.gui gui

# O ambas versiones (CLI + GUI)
make -f Makefile.gui all
```

### 3. Ejecutar
```bash
# Ejecutar GUI
make -f Makefile.gui run-gui

# O directamente
./dna_finder_gui
```

## 🚀 Uso de la Interfaz Gráfica

### Flujo Básico de Trabajo

1. **Ingreso de Secuencia**
   - Escribir directamente en el área de texto
   - O usar "Cargar archivo FASTA"

2. **Selección de Análisis**
   - **Análisis Básico**: Info básica + composición
   - **Análisis Completo**: Todos los análisis disponibles
   - **Análisis Específicos**: Traducción, ORFs, patrones

3. **Visualización de Resultados**
   - Navegar entre pestañas de resultados
   - Ver reporte completo en la última pestaña

4. **Exportación** (Opcional)
   - Usar menú "Archivo" → "Exportar resultados"
   - Guardar como archivo de texto

### Características de la Interfaz

#### Panel de Entrada
- **Validación automática**: Indica si la secuencia es válida
- **Contador de nucleótidos**: Muestra longitud en tiempo real
- **Formato automático**: Convierte a mayúsculas automáticamente

#### Panel de Análisis
- **Botones inteligentes**: Se activan solo con secuencias válidas
- **Análisis progresivo**: Puedes ejecutar análisis específicos
- **Barra de progreso**: Para análisis completos

#### Resultados
- **Pestañas organizadas**: Cada tipo de análisis en su pestaña
- **Formato legible**: Texto bien estructurado y formateado
- **Búsqueda rápida**: Fácil navegación entre resultados

## 📁 Estructura de la GUI

```
src/
├── main_gui.cpp         // Punto de entrada GUI
├── MainWindow.h         // Definición ventana principal
├── MainWindow.cpp       // Implementación GUI
├── [archivos originales]// Clases de análisis (sin cambios)
```

## 🎨 Capturas de Pantalla

La interfaz incluye:
- **Diseño moderno**: Interfaz limpia y profesional
- **Layout responsivo**: Se adapta al tamaño de ventana
- **Organización clara**: Paneles separados para entrada y resultados
- **Mensajes informativos**: Estado y progreso siempre visible

## 🔧 Solución de Problemas

### Error: "Qt development libraries not found"
```bash
# Verificar instalación Qt
pkg-config --list-all | grep -i qt

# Si no aparece nada, reinstalar Qt
# Ubuntu: sudo apt install qtbase5-dev
# macOS: brew install qt
```

### Error de compilación con MOC
```bash
# Verificar que moc esté en el PATH
which moc

# Si no está, instalar herramientas de desarrollo Qt
# Ubuntu: sudo apt install qtbase5-dev-tools
```

### Error: "Package Qt5Core not found"
```bash
# Verificar PKG_CONFIG_PATH
echo $PKG_CONFIG_PATH

# En macOS con Homebrew, puede necesitar:
export PKG_CONFIG_PATH="/opt/homebrew/lib/pkgconfig:$PKG_CONFIG_PATH"
```

## 📊 Ventajas de la Versión GUI

### Para Estudiantes
- **Visualización clara**: Resultados organizados en pestañas
- **Interacción intuitiva**: No necesita conocer comandos
- **Validación inmediata**: Retroalimentación instantánea
- **Exportación fácil**: Guardar resultados con un clic

### Para Investigadores
- **Análisis rápido**: Interfaz eficiente para análisis rutinarios
- **Carga masiva**: Soporte para archivos FASTA múltiples
- **Reportes formateados**: Resultados listos para documentación
- **Workflow visual**: Proceso de análisis claro y ordenado

## 🔄 Comparación CLI vs GUI

| Característica | CLI | GUI |
|----------------|-----|-----|
| Interfaz | Línea de comandos | Gráfica intuitiva |
| Curva de aprendizaje | Media | Baja |
| Automatización | Excelente | Limitada |
| Visualización | Texto plano | Organizada por pestañas |
| Exportación | Archivo simple | Diálogo de guardado |
| Análisis masivo | Scripts | Manual |
| Portabilidad | Alta | Requiere Qt |

## 📝 Compilación Avanzada

### Configuración de Qt Custom
```bash
# Si Qt está en ubicación no estándar
export QT_DIR=/ruta/a/qt
export PKG_CONFIG_PATH="$QT_DIR/lib/pkgconfig:$PKG_CONFIG_PATH"
```

### Compilación Debug
```bash
# Compilar con información de debug
make -f Makefile.gui gui CXXFLAGS="-std=c++11 -Wall -Wextra -g -O0 $(QT_CXXFLAGS)"
```

### Distribución
```bash
# Crear paquete distribuible
make -f Makefile.gui package
```

## 🚀 Próximas Funcionalidades

La interfaz gráfica puede expandirse con:
- 📊 Gráficos de composición nucleotídica
- 🔍 Vista de secuencia con resaltado de ORFs
- 📈 Visualización de patrones encontrados
- 🧬 Representación gráfica de la traducción
- 💾 Historial de análisis recientes
- ⚙️ Configuración de parámetros de análisis

---

¡La versión GUI hace que el análisis de ADN sea más accesible para estudiantes y investigadores!