import os
import shutil
import subprocess
import sys

# List of required tools/software
required_tools = [
    "centrifuge",
    "kraken2",
    "seqtk",
    "bwa",
    "samtools",
    "bedtools",
    "R",
    "python3",
    "git"
]

# List of Python packages with minimum versions
required_python_packages = {
    "snakemake": "7.32.4",
    "pysam": "0.6",
    "cython": None,
    "numpy": "1.24.4",
    "tqdm": None,  # No version check for tqdm
    "pandas": "2.1.2",
    "pyfaidx": "0.8.1.1",
    "biopython": "1.84",
    "scipy": "1.11.3"
}

# List of R libraries
required_r_libraries = [
    "inline",
    "gam",
    "Rcpp",
    "RcppGSL",
    "ggplot2"
]

# Map pip package names -> import names
IMPORT_NAME = {
    "biopython": "Bio",
    # most others import by the same name:
    # "snakemake": "snakemake",
    # "pysam": "pysam",
    # "numpy": "numpy",
    # "tqdm": "tqdm",
    # "pandas": "pandas",
    # "pyfaidx": "pyfaidx",
    # "scipy": "scipy",
    # "cython": "cython"
}

def check_python_package_installed(package, version=None):
    """Check if a Python package is installed and meets the version requirement."""
    modname = IMPORT_NAME.get(package, package)
    try:
        pkg = __import__(modname)
    except ImportError:
        print(f"[ERROR] Python package {package} (import '{modname}') is not installed.")
        return False

    # Try to get a version string
    installed_version = getattr(pkg, "__version__", None)
    if version and installed_version:
        try:
            # Prefer packaging, fallback to distutils if unavailable
            try:
                from packaging import version as pkg_version
            except Exception:
                from distutils.version import LooseVersion as pkg_version  # type: ignore
            if pkg_version.parse(installed_version) < pkg_version.parse(version):
                print(f"[WARNING] {package} is installed (version {installed_version}), "
                      f"but {version} or higher is required.")
                return False
        except Exception:
            # If version parsing fails, just acknowledge presence
            pass

    print(f"[OK] {package} is installed" + (f" (version {installed_version})" if installed_version else "") + ".")
    return True

def check_tool_installed(tool):
    """Check if a tool is installed and in PATH."""
    return shutil.which(tool) is not None

def check_r_library_installed(library):
    """Check if an R library is installed."""
    try:
        subprocess.run(
            ["R", "-e", f"if (!requireNamespace('{library}', quietly=TRUE)) stop('R library {library} not installed')"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )
        print(f"[OK] R library {library} is installed.")
        return True
    except subprocess.CalledProcessError:
        print(f"[ERROR] R library {library} is not installed.")
        return False

# Check tools/software
print("Checking tools and software...")
for tool in required_tools:
    if check_tool_installed(tool):
        print(f"[OK] {tool} is installed and in PATH.")
    else:
        print(f"[ERROR] {tool} is not installed or not in PATH.")

# Check Python packages
print("\nChecking Python packages...")
for package, version in required_python_packages.items():
    if not check_python_package_installed(package, version):
        continue

# Check R libraries
print("\nChecking R libraries...")
for library in required_r_libraries:
    check_r_library_installed(library)

print("\nDependency check complete!")