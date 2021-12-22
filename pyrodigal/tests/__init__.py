from . import (
    test_connection_scorer,
    test_gene,
    test_genes,
    test_nodes,
    test_pyrodigal,
    test_sequence,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_connection_scorer))
    suite.addTests(loader.loadTestsFromModule(test_gene))
    suite.addTests(loader.loadTestsFromModule(test_genes))
    suite.addTests(loader.loadTestsFromModule(test_nodes))
    suite.addTests(loader.loadTestsFromModule(test_pyrodigal))
    suite.addTests(loader.loadTestsFromModule(test_sequence))
    return suite
