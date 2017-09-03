


def latex_table(t, headers=None):
    cols = len(t)
    rows = len(t[0])
    print('\t\\begin{center}')
    print('\t\t\\begin{tabular}{' + '|c'*cols + '|}')
    print('\t\t\t\\hline')
    if headers != None:
        assert len(headers) == cols
        header_str = ''
        for h in headers:
            header_str += str(h) + '&'
        header_str = header_str[:-1] #remove trailing &
        print('\t\t\t' + header_str + '\\\\ \\hline')
    for i in range(rows):
        col = ''
        for j in range(cols):
            col += str(t[j][i]) + '&'
        col = col[:-1] #remove trailing &
        print('\t\t\t' + col + '\\\\ \\hline')
    print('\t\t\\end{tabular}')
    print('\t\\end{center}')
