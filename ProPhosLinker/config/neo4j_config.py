#!/usr/bin/env python3
"""
Neo4jConfig: Neo4j Database Connection Configuration

Configuration class for Neo4j graph database connections used in knowledge graph 
construction and functional network analysis. Provides default connection parameters 
for local Neo4j instances and supports dynamic configuration updates.

Key Features:
- Bolt protocol connection (bolt://localhost:7687)
- Default credentials (neo4j/neo4j) for development
- Optional neo4j_path for Docker/containerized deployments
- Dictionary export for easy integration with GraphDatabase.driver()
- Safe parameter updates with None handling

Usage:
    config = Neo4jConfig.to_dict()
    driver = GraphDatabase.driver(**config)
    Neo4jConfig.update_config(uri="bolt://prod-server:7687", password="secure")
"""


from dataclasses import dataclass
from typing import ClassVar


@dataclass
class Neo4jConfig:
    """Neo4j database configuration class"""
    
    uri: ClassVar[str] = "bolt://localhost:7687"
    username: ClassVar[str] = "neo4j"
    password: ClassVar[str] = "neo4j"
    neo4j_path: ClassVar[str] = ""
    
    @classmethod
    def to_dict(cls) -> dict:
        """Convert to dictionary"""
        return {
            'uri': cls.uri,
            'username': cls.username,
            'password': cls.password,
            'neo4j_path': cls.neo4j_path,
        }
    
    @classmethod
    def update_config(cls, uri: str = None, username: str = None, password: str = None, neo4j_path: str = None):
        """Update Neo4j configuration"""
        if uri:
            cls.uri = uri
        if username:
            cls.username = username
        if password:
            cls.password = password
        if neo4j_path:
            cls.neo4j_path = neo4j_path
