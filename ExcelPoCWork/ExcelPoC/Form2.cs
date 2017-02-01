using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace ExcelPoC
{
    public partial class Form2 : Form
    {
        public Form2(string[] theTables)
        {
            InitializeComponent();
            myTables = theTables;
        }

        string[] myTables;
        string myTable = string.Empty;

        public void form2_Load(object sender, EventArgs e) 
        {
            comboBox1.Items.AddRange(myTables);
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            myTable = comboBox1.SelectedItem.ToString();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            if (myTable != string.Empty)
            {
                Form1.mySelectedTable = myTable;
                this.Close();
            }
            else
            {
                MessageBox.Show("Select a Table");
            }
        }
    }
}
