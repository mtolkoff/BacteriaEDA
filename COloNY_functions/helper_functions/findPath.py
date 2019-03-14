#example path: categories['children][0]['children'][1]['children'][1]['children'][13]
#[1] returns the name of the ortholog
def findPath(KO_num, path, node):
    if(path != ''):
        token = ', '
    else:
        token = ''
    if KO_num in node['name']:
        path = list(map(int, path.split(",")))
        return([node['name'], path])
    elif 'children' in node.keys():
        for i in range(len(node['children'])):
            nodeBelow = findPath(KO_num, path + token + str(i), node['children'][i])
            if not nodeBelow == False:
                return nodeBelow
    return False
