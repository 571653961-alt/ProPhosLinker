# connect neo4j
import requests
import subprocess
import time
import os
import sys
from pathlib import Path
from neo4j import GraphDatabase
import zipfile


class Neo4jManager:
    def __init__(self, 
                http_uri="http://localhost:7474",
                bolt_uri="bolt://localhost:7687",
                username="neo4j",
                password="neo4j",
                neo4j_path=None,
                zip_relative_path = "./dataload/protein_disease_full_subgraph.zip",
                cypher_filename = "protein_disease_full_subgraph.cypher"):
        self.http_uri = http_uri
        self.bolt_uri = bolt_uri
        self.username = username
        self.password = password
        self.neo4j_path = neo4j_path
        self.zip_relative_path = zip_relative_path
        self.cypher_filename = cypher_filename
        
    def check_status(self):
        """综合检查Neo4j状态"""
        methods = [
            self._check_http,
            self._check_bolt_port,
            self._check_with_driver
        ]
        
        results = []
        for method in methods:
            try:
                results.append(method())
            except Exception as e:
                print(f"❌ 检查方法 {method.__name__} 出错: {e}\n")
                results.append(False)
        
        return any(results)
    
    def _check_http(self):
        """检查HTTP接口"""
        try:
            response = requests.get(f"{self.http_uri}/", timeout=5)
            return response.status_code == 200
        except:
            return False
    
    def _check_bolt_port(self):
        """检查Bolt端口"""
        try:
            import socket
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(2)
            result = sock.connect_ex(("localhost", 7687))
            sock.close()
            return result == 0
        except:
            return False
    
    def _check_with_driver(self):
        """使用驱动检查"""
        try:
            from neo4j import GraphDatabase
            driver = GraphDatabase.driver(self.bolt_uri, auth=(self.username, self.password))
            with driver.session() as session:
                result = session.run("RETURN 1")
                return result.single()[0] == 1
        except:
            return False
    
    def start_neo4j(self):
        """启动Neo4j"""
        print("⚠️ 启动Neo4j数据库...\n")
        
        # 尝试不同的启动方法
        methods = [
            self._start_via_command,
            self._start_via_direct_path,
            self._start_via_service
        ]
        
        for method in methods:
            if method():
                print("✅ Neo4j启动命令执行成功\n")
                return True
        
        return False
    
    def _start_via_command(self):
        """通过系统命令启动"""
        try:
            if os.name == 'nt':  # Windows
                cmd = ["neo4j", "start"]
            else:
                cmd = ["neo4j", "start"] # Linux/Mac
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            return result.returncode == 0
        except:
            return False
    
    def _start_via_service(self):
        """通过系统服务启动"""
        try:
            # Linux/Mac
            if os.name != 'nt':
                cmd = ["sudo", "systemctl", "start", "neo4j"]
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                return result.returncode == 0
            return False
        except:
            return False
    
    def _start_via_direct_path(self):
        """通过直接路径启动（兼容传入安装根目录或 bin 目录）"""
        if not self.neo4j_path:
            return False

        try:
            base = Path(self.neo4j_path)
            # 如果用户传入的是 bin 目录就直接用它，否则在安装根下寻找 bin
            bin_path = base if base.name.lower() == "bin" else base / "bin"

            if os.name == 'nt':
                exe = bin_path / "neo4j.bat"
            else:
                exe = bin_path / "neo4j start"

            # debug: 确认可执行存在
            if not exe.exists():
                print(f"Neo4j executable not found: {exe}\n", file=sys.stderr)
                return False
            cmd = [str(exe), "console"]

            # Windows: 可选使用 CREATE_NEW_CONSOLE 让 neo4j 控制台在新窗口打开
            if os.name == 'nt':
                creationflags = subprocess.CREATE_NEW_CONSOLE
                subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, creationflags=creationflags)
            else:
                subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            return True
        except Exception as e:
            print(f"_start_via_direct_path error: {e}\n", file=sys.stderr)
            return False
    
    def ensure_running(self, max_wait=90):
        """确保Neo4j运行"""
        print("************========= Neo4j数据库管理 =========************\n")
        
        # 检查当前状态
        if self.check_status():
            print("✅ Neo4j已经在运行\n")
            return True
        
        print("⚠️ Neo4j未运行，尝试启动...\n")
        
        # 启动Neo4j
        if not self.start_neo4j():
            print("❌ 无法启动Neo4j\n")
            return False
        
        # 等待启动完成
        print("⚠️ 等待Neo4j启动...\n")
        start_time = time.time()
        while time.time() - start_time < max_wait:
            if self.check_status():
                print("✅ Neo4j启动成功\n")
                return True
            time.sleep(5)
            print(".", end="", flush=True)
        
        print("❌ Neo4j启动超时\n")
        return False


    def import_data_cypher(self):
        """
        从ZIP文件导入Cypher数据
        
        Args:
            zip_relative_path: ZIP文件相对于当前脚本的路径
            cypher_filename: ZIP内的Cypher文件名（可选，不指定则自动查找）
        """
        # 获取当前脚本所在目录
        current_dir = os.path.dirname(os.path.abspath(__file__))
        zip_path = os.path.join(current_dir, self.zip_relative_path)
        
        # 确保数据库运行
        if not self.ensure_running():
            print("❌ 无法启动Neo4j，导入中止")
            return False
        
        try:
            print(f"📦 开始导入数据从: {self.zip_relative_path}")
            
            # 解压ZIP文件
            with zipfile.ZipFile(zip_path, 'r') as z:
                # 如果未指定cypher文件名，自动查找
                if cypher_filename is None:
                    file_list = z.namelist()
                    cypher_files = [f for f in file_list if f.endswith('.cypher')]
                    if not cypher_files:
                        print("❌ 未找到.cypher文件")
                        return False
                    cypher_filename = cypher_files[0]
                
                # 解压到当前目录
                z.extract(cypher_filename, current_dir)
                extracted_path = os.path.join(current_dir, cypher_filename)
                print(f"✅ 文件解压成功: {cypher_filename}")
            
            # 执行Cypher文件
            print("🚀 执行Cypher导入...")
            cmd = f'cypher-shell -u {self.username} -p {self.password} -f "{extracted_path}"'
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                print("✅ Cypher导入成功!")
                if result.stdout.strip():
                    print("输出:", result.stdout.strip())
            else:
                print("❌ Cypher导入失败!")
                print("错误:", result.stderr)
                return False
            
            return True
            
        except Exception as e:
            print(f"❌ 导入过程中出错: {e}")
            return False
        finally:
            # 清理临时文件
            if 'extracted_path' in locals() and os.path.exists(extracted_path):
                os.remove(extracted_path)
                print(f"🧹 已清理临时文件: {cypher_filename}")
# 使用示例
if __name__ == "__main__":
    manager = Neo4jManager(
        username="neo4j",
        password="neo4j",
        neo4j_path="D:/software/neo4j-community-4.2.3/bin"  # 可选
    )
    
    if manager.ensure_running():
        print("✅ Neo4j已就绪，可以开始工作\n")
    else:
        print("❌ Neo4j启动失败，请手动检查\n")
        sys.exit(1)