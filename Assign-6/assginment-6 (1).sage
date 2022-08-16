def coppersmith_howgrave_univariate(pol, modulus, beta, mm, tt, XX):
    
    dd = pol.degree()
    polZ = pol.change_ring(ZZ)
    
    nn = dd * mm + tt    
    x = polZ.parent().gen()

   
    gg = []
    BB = Matrix(ZZ, nn)
    for ii in range(mm):
        for jj in range(dd):
            gg.append((x * XX)**jj * modulus**(mm - ii) * polZ(x * XX)**ii)
    for ii in range(tt):
        gg.append((x * XX)**ii * polZ(x * XX)**mm)
 

    for ii in range(nn):
        for jj in range(ii+1):
            BB[ii, jj] = gg[ii][jj]

    
    BB = BB.LLL()
    
    new_pol = 0
    for ii in range(nn):
        new_pol += x**ii * BB[0, ii] / XX**ii

    potential_roots = new_pol.roots()

    roots = []
    for root in potential_roots:
        if root[0].is_integer():
            result = polZ(ZZ(root[0]))
            if gcd(modulus, result) >= modulus^beta:
                roots.append(ZZ(root[0]))

    return roots

e = 5
N = 84364443735725034864402554533826279174703893439763343343863260342756678609216895093779263028809246505955647572176682669445270008816481771701417554768871285020442403001649254405058303439906229201909599348669565697534331652019516409514800265887388539283381053937433496994442146419682027649079704982600857517093
C = 16406040638929092928886376275629768360969123830963841491851182018240415264239447408750645116516593638454711155935903116677173236783008186384387499645921457504138916204718432271060637800126956587541917183968180412666152342086100869010734708863775005626622338061044022246658000252552309180133555987688387216897

# RSA known parameters
ZmodN = Zmod(N);

def break_RSA(p_str, max_length_M):
    global e, C, ZmodN
    p_binary_str = ''.join(['{0:08b}'.format(ord(x)) for x in p_str])

    for length_M in range(0, max_length_M+1, 4):        
   
        P.<M> = PolynomialRing(ZmodN)
        pol = ((int(p_binary_str, 2)<<length_M) + M)^e - C
        dd = pol.degree()

        # Tweak those
        beta = 1                                                          
        tt = floor(dd * (ceil(beta**2 / (dd * (beta/7)))) * ((1/beta) - 1))    
        XX = ceil(N**((beta**2/dd) - (beta/7)))  

        roots = coppersmith_howgrave_univariate(pol, N, beta, (ceil(beta**2 / (dd * (beta/7)))), tt, XX)

        if roots:
            return '{0:b}'.format(roots[0])

    print('No solution found\n')
    return 0

def bin2char(binary_message):
    message = ""
    if(len(binary_message)%8 != 0):
        binary_message = (8-len(binary_message)%8)*"0" + binary_message
        for i in range(0,len(binary_message),8):
            decimal_number = int(binary_message[i:i+8],2)
            message += chr(decimal_number)
    else:
        for i in range(0,len(binary_message),8):
            decimal_number = int(binary_message[i:i+8],2)
            message += chr(decimal_number)  

    return message

if __name__ == "__main__":
    #binary_message = break_RSA("Achievers: This door has RSA encryption with exponent 5 and the password is ", 300 )
    if(break_RSA("Achievers: This door has RSA encryption with exponent 5 and the password is ", 300 )!=0):
        print("The root in binary is: " + break_RSA("Achievers: This door has RSA encryption with exponent 5 and the password is ", 300 ))
        print("Password : " + bin2char(break_RSA("Achievers: This door has RSA encryption with exponent 5 and the password is ", 300 )))