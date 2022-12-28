from dnacurve import CurvedDNA
from dnacurve import Sequence
import os
test = Sequence('workflow/scripts/TSS_seq_LUZ7_enriched_t5.plus.fa', name="test")
print(test)
result = CurvedDNA(test, 'trifonov', name='Example')
print(result.curvature[:,])
result.plot("out.png")