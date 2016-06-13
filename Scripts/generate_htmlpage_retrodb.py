'''
@author: Andrea Caddedu, Liz Wylie
generates a very basic html file with a table of all the reactions in db
a typical html table is expressed:

<table border="1">
<tr>
<td>row 1, cell 1</td>
<td>row 1, cell 2</td>
</tr>
<tr>
<td>row 2, cell 1</td>
<td>row 2, cell 2</td>
</tr>
</table> 


a typical db entry:
 u'condition': u'Co2(CO)6(P(OPh)3)2.Co2(CO)8.toluene',
 u'name': u'Pauson-Khand Reaction',
 u'popularity_score': 0.7269344658439667,
 u'products': [u'O=C1C([#6])=CC([#6])C1'],
 u'reaction_smarts': u'[C:6]1[C:5]([#6:10])[C:3]=[C:2]([#6:4])[C:8]1=[O:9]>>[#6:4][C:2]#[C:3].[C:6]=[C:5][#6:10].[C-:8]#[O+:9]',
 u'synthon_count'

'''

import pymongo
from rdkit import Chem


if __name__ == '__main__':
    db = pymongo.MongoClient('localhost').new_data

    with open('reactionlist.html', "w") as f:
        # this is the start of the html header. It's required by the browser!
        f.write('<html> \n<head>')
        # this is just a css for making the output less ugly
        f.write(
            "<meta charset=""utf-8""> <title>Chematicaweb0.0</title> <link href=""http://twitter.github.com/bootstrap/1.4.0/bootstrap.css"" rel=""stylesheet""><style>.content {padding-top: 80px;}</style>")
        # this closes the html head and start the body. Required
        f.write('<br></head>\n<body>')
        f.write('<table border="1">')  # start table
        f.write(
            '<tr><td> Reaction Name</td><td> Reaction Smarts</td><td>Products</td><td>retrons img</td><td>====></td><td>synthons img</td></tr>')
        for entry in db.retro.find():
            # save two png file, and add a line in html table
            img = str(entry['_id']) + ".png"
            image_size = (200, 150)

            react, prod = entry['reaction_smarts'].encode('ascii').split('>>')
            reacts = react.split('.')
            prods = prod.split('.')

            try:
                reactants = [Chem.MolFromSmarts(a) for a in reacts]
            except:
                print "Error in parsing reactants for reaction #" + str(entry['_id'])
                print 'reactants:', reactants

            try:
                products = [Chem.MolFromSmarts(b) for b in prods]
            except:
                print "Error in parsing products for reaction #" + str(entry['_id'])
                print 'products:', products

            map(lambda m: m.UpdatePropertyCache(strict=False), reactants)
            map(lambda m: m.UpdatePropertyCache(strict=False), products)

            try:
                reactant_images = []
                for i, motif in enumerate(reactants):
                    rimg = str("r" + i + img)
                    Chem.MolToFile(
                        motif, rimg, size=image_size, imageType="png")
                    reactant_images.append(rimg)
            except:
                print "Error in creating reactant images for reaction #" + str(entry['_id'])

            try:
                product_images = []
                for i, motif in enumerate(products):
                    pimg = str("p" + i + img)
                    Chem.MolToFile(
                        motif, pimg, size=image_size, imageType="png")
                    product_images.append(pimg)
            except:
                print "Error in creating product images for reaction #" + str(entry['_id'])

            entry_name = entry['name'].encode('ascii')
            rs = entry['reaction_smarts'].encode('ascii')
            productsmarts = [i.encode('ascii') for i in entry['products']]
            productstring = productsmarts.pop()
            while productsmarts != []:
                smarts = productsmarts.pop()
                productstring = productstring + ', ' + smarts

            a = '<tr><td>' + entry_name + '</td><td>' + \
                rs + '</td><td>' + productstring
            a = a + '</td><td> <img src=' + \
                reactant_images[
                    0] + '></td><td> ===> </td><td> <img src=' + product_images[0] + '></td></tr>'

            f.write(a)

        f.write('</table></html>')
