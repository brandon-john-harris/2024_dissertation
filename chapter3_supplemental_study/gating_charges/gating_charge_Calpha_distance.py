from Bio.PDB import PDBParser


def calc_ca_distance(structure, ref_chain_id, ref_residue_id, test_residue_ids, insertion_code=''):
    distances={}

    ref_chain = structure[0][ref_chain_id]
    for residue in ref_chain:
        if residue.id[1] == ref_residue_id:
            ref_residue = residue
            break

    if ref_residue is None:
        raise ValueError("reference residue not found")

    ref_ca = ref_residue['CA']

    for test_residue_id in test_residue_ids:
        test_residue = None
        for residue in ref_chain:
            if residue.id[1] == test_residue_id:
                test_residue = residue
                break

        if test_residue is None:
            raise ValueError(f'test residue {test_residue_id} not found in structure')

        test_ca = test_residue['CA']
        distance = ref_ca - test_ca
        distances[test_residue_id] = distance

    return distances


# Load structure
root_path = '/Users/user/Data/'
input_data = ['6uz0.pdb']  # adjust as needed
for input_file in input_data:
    Chai1_model = f'{root_path}/6uz0/reference/{input_file}'  # adjust as needed
    parser = PDBParser()
    print('----------------------')
    print(f'File: {input_file}')
    structure = parser.get_structure('Cha1_model', Chai1_model)

    ref_chain_id = 'A'  # adjust as needed
    # adjust as needed: HCS: [gating charges]
    residue_coordinates = {'VSD1': (169, [220, 223, 226, 229]),
                           'VSD2': (761, [809, 812, 815, 818, 821]),
                           'VSD3': (1252, [1302, 1305, 1308, 1311, 1314, 1318]),
                           'VSD4': (1573, [1625, 1628, 1631, 1634, 1637, 1640])}
    for vsd, (ref_HCS_residue_id, gating_charge_residue_ids) in residue_coordinates.items():
        print(vsd)
        distances = calc_ca_distance(structure, ref_chain_id, ref_HCS_residue_id, gating_charge_residue_ids, ' ')
        for resid, distance in distances.items():
            print(distance)

exit()
