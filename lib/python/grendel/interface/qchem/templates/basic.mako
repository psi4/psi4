<%
    # shorter names for some stuff...
    mol = computation.molecule
    det = computation.details
%>
$molecule
  ${mol.charge} ${mol.multiplicity}
  ${mol.xyz_string(header=False)}
$end

<%
    optional = [
        'mem_static',
        'mem_total',
        'scf_algorithm',
        'xc_grid',
        'max_scf_cycles'
    ]
%>

$rem
% if det.method == Methods.DFT:
  exchange = ${det.exchange}
% elif det.method == Methods.HF:
  exchange = HF
% endif
  basis = ${basis}
% if computation.needs_property('EnergyGradient'):
  jobtype = force
% elif computation.needs_property('Hessian'):
  jobtype = freq
% else:
  jobtype = sp
% endif
  ${det.keywordify('{key} = {value}', *optional)}
$end
