from vamdclib import nodes
from vamdclib import request as r
from vamdclib import specmodel as m


nl = nodes.Nodelist()
nl.findnode('cdms')
cdms = nl.findnode('cdms')

request = r.Request(node=cdms)


# Retrieve all species from CDMS
result = request.getspecies()
molecules = result.data['Molecules']

h2co = [x for x in molecules.values()
        if (hasattr(x,'MolecularWeight') and
            (x.StoichiometricFormula)==('CH2O') and
            (x.MolecularWeight=='30'))][0]

h2co_inchikey = h2co.InChIKey

# query everything for h2co
query_string = "SELECT ALL WHERE VAMDCSpeciesID='%s'" % h2co.VAMDCSpeciesID
request.setquery(query_string)
result = request.dorequest()

def Qrot(temperature):
    Q = m.calculate_partitionfunction(result.data['States'],
                                      temperature=temperature)[h2co.Id]
    return Q
