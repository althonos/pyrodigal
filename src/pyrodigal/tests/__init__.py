from . import (
    test_connection_scorer,
    test_gene,
    test_genes,
    test_gene_finder,
    test_mask,
    test_nodes,
    test_sequence,
    test_training_info,
    test_doctest,
)


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_connection_scorer))
    suite.addTests(loader.loadTestsFromModule(test_gene))
    suite.addTests(loader.loadTestsFromModule(test_genes))
    suite.addTests(loader.loadTestsFromModule(test_gene_finder))
    suite.addTests(loader.loadTestsFromModule(test_mask))
    suite.addTests(loader.loadTestsFromModule(test_nodes))
    suite.addTests(loader.loadTestsFromModule(test_sequence))
    suite.addTests(loader.loadTestsFromModule(test_training_info))
    test_doctest.load_tests(loader, suite, pattern)
    return suite
