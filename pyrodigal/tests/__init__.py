from . import test_gene, test_genes, test_nodes, test_pyrodigal, test_connection_scorer

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_gene))
    suite.addTests(loader.loadTestsFromModule(test_genes))
    suite.addTests(loader.loadTestsFromModule(test_nodes))
    suite.addTests(loader.loadTestsFromModule(test_pyrodigal))
    suite.addTests(loader.loadTestsFromModule(test_connection_scorer))
    return suite
