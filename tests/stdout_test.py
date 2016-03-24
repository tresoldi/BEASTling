from __future__ import unicode_literals
import sys
from contextlib import contextmanager

from six import BytesIO

import beastling.beastxml
import beastling.configuration
from .util import WithConfigAndTempDir


@contextmanager
def capture(func, *args, **kw):
    out, sys.stdout = sys.stdout, BytesIO()
    func(*args, **kw)
    sys.stdout.seek(0)
    yield sys.stdout.read()
    sys.stdout = out


class Tests(WithConfigAndTempDir):
    def test_stdout(self):
        xml = beastling.beastxml.BeastXml(self.make_cfg("tests/configs/basic.conf"))
        with capture(xml.write_file, 'stdout') as output:
            self.assertIn('<?xml version=', output.decode('utf8'))
