#!/usr/bin/env python3
"""
ProPhosLinker Setup Script with Custom Installation

Comprehensive setup script for ProPhosLinker - Protein-Phosphoprotein network analysis 
pipeline. Features automated environment validation (Python/R/Neo4j), intelligent R package 
installation (37 packages), and optional Neo4j data import.

Key Features:
- Python 3.9+ environment validation
- R installation detection and 37-package auto-installation  
- Neo4j binary detection (neo4j/cypher-shell commands)
- Cross-platform Windows/Linux/macOS compatibility
- Production-ready dependency management

Installation Flow:
1. Environment validation (Python/R/Neo4j)  
2. Python package installation via pip
3. R package installation (skip if exists)
4. Optional Neo4j database population
"""

from setuptools import setup, find_packages
import os
from setuptools.command.install import install
import subprocess
import sys
import zipfile
import getpass


class CustomInstallCommand(install):
    """Custom installation command with environment detection and Neo4j data import"""
    
    def run(self):
        """Main installation workflow"""
        print("\n" + "="*60)
        print("🔍 Starting environment detection and installation...")
        
        # 1. Check Python version
        if not self.check_python_version():
            sys.exit(1)
        
        # 2. Check R installation
        if not self.check_r_installation():
            sys.exit(1)
        
        # 3. Check Neo4j installation  
        if not self.check_neo4j_installation():
            print("⚠️  Neo4j not installed or version too old, please install Neo4j 4.0+ manually")
            # Continue installation as user may install Neo4j later
        
        # 4. Execute standard Python package installation
        print("\n📦 Installing Python dependencies...")
        install.run(self)
        
        # 5. Install R packages
        print("\n📊 Installing R dependencies...")
        self.install_r_packages()
        
        print("\n🎉 Installation completed!")

    def check_python_version(self):
        """Check Python version compatibility"""
        print("\n1. Checking Python version...")
        required_version = (3, 9)
        current_version = sys.version_info[:2]
        
        if current_version < required_version:
            print(f"❌ Python version too low: {sys.version}")
            print(f"   Requires Python {required_version[0]}.{required_version[1]}+")
            return False
        else:
            print(f"✅ Python version OK: {sys.version}")
            return True

    def check_r_installation(self):
        """Check R installation availability"""
        print("\n2. Checking R installation...")
        try:
            # Try running R command
            result = subprocess.run(
                ["R", "--version"], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            
            if result.returncode == 0:
                # Extract version information
                version_line = result.stdout.split('\n')[0]
                print(f"✅ R installed successfully: {version_line}")
                return True
            else:
                print("❌ R command execution failed")
                return False
                
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
            print("❌ R not found, please install R first: https://www.r-project.org/")
            return False

    def check_neo4j_installation(self):
        """Check Neo4j installation status"""
        print("\n3. Checking Neo4j installation...")
        
        # Method 1: Check neo4j command
        try:
            result = subprocess.run(
                ["neo4j", "version"], 
                capture_output=True, 
                text=True, 
                timeout=10,
                shell=True  # Required for Windows compatibility
            )
            if result.returncode == 0:
                version_line = result.stdout.split('\n')[0]
                print(f"✅ Neo4j installed: {version_line}")
                return True
        except Exception as e:
            print(f"neo4j command failed: {e}")
        
        # Method 2: Check cypher-shell command
        try:
            result = subprocess.run(
                "cypher-shell --version",  # String command for shell compatibility
                capture_output=True, 
                text=True, 
                timeout=10,
                shell=True  # Essential for Windows
            )
            if result.returncode == 0:
                print(f"✅ Cypher-shell available: {result.stdout.strip()}")
                return True 
            else:
                print(f"cypher-shell non-zero exit: {result.stderr}")
        except Exception as e:
            print(f"cypher-shell command error: {e}")
        
        print("❌ Neo4j installation not detected")
        return False

    def install_r_packages(self):
        """Install required R packages"""
        r_packages = [
            "optparse",
            "readr",
            "stringr",
            "dplyr",
            "vegan",
            "ggrepel",
            "ggplot2",
            "NMF",
            "tidyverse",
            "doParallel",
            "WGCNA",
            "patchwork",
            "pheatmap",
            "plyr",
            "viridis",
            "grid",
            "flashClust",
            "ggsankeyfier",
            "limma",
            "statmod",
            "colorspace",
            "Mfuzz",
            "reshape2",
            "tibble",
            "igraph",
            "ggraph",
            "tidygraph",
            "tidyr",
            "ggforce",
            "ggpubr",
            "clusterProfiler",
            "org.Hs.eg.db",
            "enrichplot",
            "hmisc",
            "bootnet",
            "graphlayouts",
            "scatterpie",
            "ggsci",
            "ggnewscale",
            "svglite",
            "ggiraph"
        ]
        
        print(f"📥 Installing {len(r_packages)} R packages...")
        
        installed_count = 0
        for package in r_packages:
            try:
                # Check if package installed and install if missing
                install_cmd = f'R -e "if(!require(\'{package}\', quietly=TRUE)) {{ install.packages(\'{package}\', repos=\'https://cloud.r-project.org/\') }}; quit(status=0)"'
                result = subprocess.run(install_cmd, shell=True, capture_output=True, text=True, timeout=120)

                if result.returncode == 0:
                    print(f"   ✅ {package} ready")
                    installed_count += 1
                else:
                    # Fallback installation
                    install_cmd = f'R -e "install.packages(\'{package}\', repos=\'https://cloud.r-project.org/\')"'
                    install_result = subprocess.run(
                        install_cmd, shell=True, capture_output=True, text=True, timeout=120
                    )
                    
                    if install_result.returncode == 0:
                        print(f"   ✅ {package} installed successfully")
                        installed_count += 1
                    else:
                        print(f"   ❌ {package} installation failed")
                        
            except Exception as e:
                print(f"   ⚠️  {package} installation error: {e}")
        
        print(f"📊 R package installation complete: {installed_count}/{len(r_packages)}")

    def check_neo4j_service(self):
        """Check if Neo4j service is running"""
        try:
            # Check Bolt port
            import socket
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(5)
            result = sock.connect_ex(("localhost", 7687))
            sock.close()
            
            if result == 0:
                print("✅ Neo4j service running")
                return True
            else:
                print("❌ Neo4j service not running")
                return False
                
        except Exception as e:
            print(f"❌ Neo4j service check error: {e}")
            return False


# Read README file
try:
    with open('README.md', 'r', encoding='utf-8') as f:
        long_description = f.read()
except:
    long_description = "A package for Protein-Phosphorylation Site Correlation Analysis"

# Read requirements.txt
try:
    with open('requirements.txt', 'r', encoding='utf-8') as f:
        requirements = f.read().splitlines()
except:
    requirements = [
        "pandas>=2.0.0",
        "numpy>=1.24.0", 
        "scipy>=1.10.0",
        "scikit-learn>=1.3.0",
        "statsmodels>=0.14.0",
        "requests>= 2.25.0",
        
        # Visualization
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
        
        # Network Analysis
        "networkx>=3.1",
        "python-igraph>=0.11.0",
        "leidenalg>=0.10.0",
        "igraph>=0.11.0",
        
        # Database
        "neo4j>=5.15.0",
        
        # Configuration & Utilities
        "pyyaml>=6.0",
        "jinja2>=3.1.0",
        
        # R Integration
        "rpy2>=3.5.10",
        
        # Web Automation
        "playwright>=1.40.0",
        
        # Development
        "setuptools>=65.0.0"
    ]

setup(
    name="ProPhosLinker",
    version="0.1.0",
    author="Pan Ying",
    author_email="panying2@genomics.cn",
    description="A package for Protein - Phosphoprotein network analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/ProPhosLinker",
    packages=find_packages(),
    package_data={
        'ProPhosLinker': [
            'database/dataload/*.zip',
            'scripts/*.R',
            'scripts/**/*',
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "ProPhosLinker=ProPhosLinker.cli:main"
        ],
    },
    include_package_data=True,
    setup_requires=['setuptools', 'wheel'],
    cmdclass={
        'install': CustomInstallCommand,
    },
)
