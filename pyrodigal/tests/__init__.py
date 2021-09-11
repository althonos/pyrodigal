from . import test_gene, test_genes, test_nodes, test_pyrodigal

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_gene))
    suite.addTests(loader.loadTestsFromModule(test_genes))
    suite.addTests(loader.loadTestsFromModule(test_nodes))
    suite.addTests(loader.loadTestsFromModule(test_pyrodigal))
    return suite
