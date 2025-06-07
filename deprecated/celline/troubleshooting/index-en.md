# Troubleshooting

This document explains problems that may occur while using Celline and their solutions.
It's a comprehensive troubleshooting guide covering common issues to advanced debugging methods.

## 🚨 Common Issues and Solutions

### Installation Related

#### Issue: pip install fails

**Symptoms:**
```bash
pip install celline
ERROR: Could not find a version that satisfies the requirement celline
```

**Solutions:**

```bash
# 1. Update pip and setuptools to latest
pip install --upgrade pip setuptools

# 2. Check Python version
python --version  # 3.10+ required

# 3. Use virtual environment
python -m venv celline-env
source celline-env/bin/activate  # Linux/Mac
pip install celline

# 4. Install development version
pip install git+https://github.com/YUYA556223/celline.git
```

#### Issue: Dependency conflicts

**Symptoms:**
```bash
ERROR: pip's dependency resolver does not currently have a strict priority resolver
```

**Solutions:**

```bash
# 1. Reinstall in clean environment
pip uninstall celline
pip cache purge
pip install celline

# 2. Specify particular version
pip install 'celline==0.1.10'

# 3. Use --force-reinstall
pip install --force-reinstall celline
```

### Configuration Related

#### Issue: R environment not found

**Symptoms:**
```bash
celline run create_seurat
Error: R path not found
```

**Solutions:**

```bash
# 1. Check if R is installed
which R
R --version

# 2. Manually set R path
celline config --r-path /usr/bin/R

# 3. Set environment variables
export R_HOME=/usr/lib/R
export PATH=$PATH:/usr/bin

# 4. Manual editing of setting.toml
[R]
r_path = "/opt/R/4.3.0/bin/R"
```

#### Issue: Configuration file cannot be loaded

**Symptoms:**
```bash
celline: command not found
# or
Error: setting.toml not found
```

**Solutions:**

```bash
# 1. Check if in correct directory
pwd
ls -la setting.toml

# 2. Initialize project
celline init

# 3. Check configuration file permissions
ls -la setting.toml
chmod 644 setting.toml

# 4. Regenerate configuration
rm setting.toml
celline init
```

### Data Acquisition Related

#### Issue: Sample addition fails

**Symptoms:**
```bash
celline run add GSE123456
Error: No resolver found for GSE123456
```

**Solutions:**

```bash
# 1. Check network connection
ping www.ncbi.nlm.nih.gov

# 2. Synchronize database
celline run sync_DB

# 3. Verify sample ID
# GSE: Gene Expression Omnibus Study
# GSM: Gene Expression Omnibus Sample  
# SRR: Sequence Read Archive Run

# 4. Manual sample verification
curl "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123456"

# 5. Execute in debug mode
celline --debug run add GSE123456
```

#### Issue: Download interrupted

**Symptoms:**
```bash
celline run download
Download failed: Connection timeout
```

**Solutions:**

```bash
# 1. Check network settings
curl -I https://ftp.ncbi.nlm.nih.gov/

# 2. Adjust timeout settings
# Increase wait_time in setting.toml
[fetch]
wait_time = 10

# 3. Reduce parallelism
celline run download --nthread 1

# 4. Proxy settings (if needed)
export http_proxy=http://proxy.example.com:8080
export https_proxy=http://proxy.example.com:8080

# 5. Partial download
# Re-execute only failed samples
```

### Processing Related

#### Issue: Cell Ranger not found

**Symptoms:**
```bash
celline run count
Error: cellranger: command not found
```

**Solutions:**

```bash
# 1. Check Cell Ranger installation
which cellranger
cellranger --version

# 2. Add to path
export PATH=$PATH:/opt/cellranger-7.0.0

# 3. Persistent path setting
echo 'export PATH=$PATH:/opt/cellranger-7.0.0' >> ~/.bashrc
source ~/.bashrc

# 4. Download Cell Ranger
# https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

# 5. Use alternative methods
# Consider processing methods without Cell Ranger
```

#### Issue: Out of memory error

**Symptoms:**
```bash
celline run preprocess
Error: MemoryError: Unable to allocate array
```

**Solutions:**

```bash
# 1. Check available memory
free -h
top

# 2. Reduce parallelism
celline config --nthread 1

# 3. Add swap file (Linux)
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

# 4. Reduce sample size
# Split large datasets for processing

# 5. Change processing method
# Use memory-efficient processing options
```

### Execution Environment Related

#### Issue: Parallel processing not working

**Symptoms:**
```bash
celline run download --nthread 4
# Only processed one at a time
```

**Solutions:**

```bash
# 1. Check execution system
celline config
# Confirm system = "multithreading"

# 2. Check available CPUs
nproc
lscpu

# 3. Update configuration
celline config --system multithreading --nthread 4

# 4. Check Python environment
python -c "import concurrent.futures; print('OK')"

# 5. Check with debug output
celline --verbose run download --nthread 4
```

#### Issue: PBS/Slurm jobs not submitted

**Symptoms:**
```bash
celline run count
Error: PBS command not found
```

**Solutions:**

```bash
# 1. Check cluster environment
qstat  # PBS
squeue  # Slurm

# 2. Load modules
module load pbs
module load slurm

# 3. Check configuration
celline config --system PBS --pbs-server cluster-name

# 4. Check job templates
ls templates/controllers/
cat templates/controllers/PBS.sh

# 5. Test manual job submission
qsub -l select=1:ncpus=1 -l walltime=00:10:00 << EOF
echo "Test job"
EOF
```

## 🔍 Debug Mode

### Enable Detailed Logging

```bash
# 1. Execute with debug flag
celline --debug run preprocess

# 2. Check verbose logs
celline --verbose run download

# 3. Check log files
tail -f celline.log

# 4. Set log level for specific modules
export CELLINE_LOG_LEVEL=DEBUG
export CELLINE_LOG_FILE=debug.log
```

### Log File Locations

```bash
# Common log file locations
./celline.log                    # Current directory
~/.celline/logs/                 # Home directory
resources/*/log/                 # Sample-specific logs
/tmp/celline_*.log              # Temporary log files
```

### Using Python Debugger

```python
# Python script for debugging
import pdb
from celline import Project
from celline.functions.preprocess import Preprocess

# Start debugger
pdb.set_trace()

project = Project("./")
preprocess = Preprocess()
result = preprocess.call(project)
```

## 📊 Performance Issues

### Memory Usage Monitoring

```bash
# 1. Monitor system resources
htop
top
free -h

# 2. Memory usage by process
ps aux | grep celline
pmap $(pgrep -f celline)

# 3. Python memory profiler
pip install memory-profiler
python -m memory_profiler script.py
```

### Processing Time Measurement

```bash
# 1. Execute with time measurement
time celline run preprocess

# 2. Detailed profiling
python -m cProfile -o profile.stats script.py
python -c "import pstats; pstats.Stats('profile.stats').sort_stats('time').print_stats(20)"

# 3. Time measurement for each step
celline --time run preprocess
```

### Disk Usage Check

```bash
# 1. Check disk usage
df -h
du -sh resources/ data/ results/

# 2. Identify large files
find . -type f -size +1G

# 3. Clean up unnecessary files
# Delete temporary files
find resources/*/tmp -type f -delete

# Rotate log files
find . -name "*.log" -mtime +7 -delete
```

## 🔧 Configuration Issues

### Configuration File Validation

```bash
# 1. Check setting.toml syntax
python -c "import toml; toml.load('setting.toml')"

# 2. Check configuration values
celline config --show

# 3. Restore default configuration
mv setting.toml setting.toml.backup
celline init

# 4. Gradual configuration restoration
# Create minimal setting.toml with basic settings only
```

### Environment Variable Issues

```bash
# 1. Check environment variables
env | grep -i celline
env | grep -i python
env | grep -i r_

# 2. Check paths
echo $PATH
echo $PYTHONPATH
echo $R_HOME

# 3. Reset environment variables
unset CELLINE_CONFIG
unset R_HOME
```

## 🌐 Network Issues

### Connection Problem Diagnosis

```bash
# 1. Basic connection test
ping www.ncbi.nlm.nih.gov
curl -I https://ftp.ncbi.nlm.nih.gov/

# 2. Check DNS resolution
nslookup www.ncbi.nlm.nih.gov
dig www.ncbi.nlm.nih.gov

# 3. Check proxy settings
echo $http_proxy
echo $https_proxy

# 4. Check firewall settings
# Administrator consultation required
```

### SSL/TLS Certificate Issues

```bash
# 1. Check certificates
openssl s_client -connect www.ncbi.nlm.nih.gov:443

# 2. Ignore certificates (temporary)
export PYTHONHTTPSVERIFY=0
pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org celline

# 3. Update certificates
# Update system certificates
```

## 🔄 Data Integrity Issues

### Database Repair

```bash
# 1. Check database files
ls -la src/celline/DB/*.parquet

# 2. Re-synchronize database
celline run sync_DB

# 3. Delete corrupted files
rm src/celline/DB/*.parquet
celline run sync_DB

# 4. Manual database verification
python -c "
import polars as pl
df = pl.read_parquet('src/celline/DB/SRA_GSM.parquet')
print(f'Records: {len(df)}')
print(df.head())
"
```

### File System Issues

```bash
# 1. Check file system
fsck /dev/sda1  # Requires administrator privileges

# 2. Check disk space
df -h
du -sh ./

# 3. Check inode usage
df -i

# 4. Fix file permissions
chmod -R 755 resources/
chmod -R 644 *.toml *.log
```

## 🆘 Emergency Procedures

### Safe Process Termination

```bash
# 1. Stop with Ctrl+C
# Press Ctrl+C in terminal

# 2. Force kill process
ps aux | grep celline
kill -9 <PID>

# 3. Stop background jobs
jobs
kill %1  # Stop by job number

# 4. Reduce system load
killall python
killall R
```

### Data Recovery

```bash
# 1. Restore from backup
cp -r backup/data/ ./
cp backup/setting.toml ./

# 2. Partial recovery
# Delete only failed samples
rm -rf resources/GSM123456/
rm -rf data/GSM123456/

# 3. Gradual recovery
# Process samples one by one
for sample in GSM123456 GSM123457; do
    celline run preprocess --sample $sample
done
```

### Log File Inspection Methods

```bash
# 1. Search for error messages
grep -i error celline.log
grep -i failed resources/*/log/*.log

# 2. Check logs chronologically
tail -f celline.log

# 3. Check logs for specific period
grep "2024-01-15" celline.log

# 4. Compress and save log files
gzip old_logs/*.log
tar -czf logs_backup.tar.gz logs/
```

## 📞 Support Resources

### Information Gathering Commands

```bash
#!/bin/bash
# debug_info.sh - Collect debug information

echo "=== Celline Debug Information ==="
echo "Date: $(date)"
echo "OS: $(uname -a)"
echo "Python: $(python --version)"
echo

echo "=== Celline Version ==="
python -c "import celline; print(f'Celline: {celline.__version__}')" 2>/dev/null || echo "Celline not installed"
echo

echo "=== Environment ==="
echo "PATH: $PATH"
echo "PYTHONPATH: $PYTHONPATH"
echo "R_HOME: $R_HOME"
echo

echo "=== System Resources ==="
echo "Memory:"
free -h
echo "Disk:"
df -h
echo "CPU:"
nproc
echo

echo "=== Configuration ==="
if [ -f setting.toml ]; then
    echo "setting.toml found"
    cat setting.toml
else
    echo "setting.toml not found"
fi

if [ -f samples.toml ]; then
    echo "samples.toml found"
    wc -l samples.toml
else
    echo "samples.toml not found"
fi
echo

echo "=== Recent Errors ==="
if [ -f celline.log ]; then
    echo "Last 10 error lines:"
    grep -i error celline.log | tail -10
else
    echo "No log file found"
fi
```

### Community Resources

- **GitHub Issues**: https://github.com/YUYA556223/celline/issues
- **Discussions**: https://github.com/YUYA556223/celline/discussions
- **Documentation**: This documentation
- **Example Projects**: https://github.com/YUYA556223/celline-examples

### Creating Bug Reports

When reporting issues, please include the following information:

```markdown
## Environment
- OS: (e.g., Ubuntu 20.04, macOS 13.0, Windows 11)
- Python Version: (e.g., 3.11.0)
- Celline Version: (e.g., 0.1.10)
- Installation Method: (e.g., pip, uv, from source)

## Expected Behavior
Describe the expected behavior

## Actual Behavior
Describe the actual behavior

## Steps to Reproduce
1. 
2. 
3. 

## Error Messages
```
Paste error messages here
```

## Additional Context
- Configuration file (setting.toml)
- Log files (relevant parts only)
- Other relevant information
```

---

> **Info**: If these solutions don't resolve your issue, please report the problem with detailed information at [GitHub Issues](https://github.com/YUYA556223/celline/issues).