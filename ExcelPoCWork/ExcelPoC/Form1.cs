using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Data.Odbc;

namespace ExcelPoC
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        public static string mySelectedTable = string.Empty;

        private void button2_Click(object sender, EventArgs e)
        {
            OpenFileDialog myOpenFileDiag = new OpenFileDialog();
            myOpenFileDiag.Title = "Select Build Sheet File";
            //TODO: alter to being current user's documents folder
            myOpenFileDiag.InitialDirectory = @"C:\Users";
            myOpenFileDiag.FileName = textBox1.Text;
            myOpenFileDiag.Filter = "Excel 97-03 Sheet(*.xls)|*.xls|Excel Sheet(*.xlsx)|*.xlsx|All Files(*.*)|*.*";
            myOpenFileDiag.FilterIndex = 2;
            myOpenFileDiag.RestoreDirectory = true;
            if (myOpenFileDiag.ShowDialog() == DialogResult.OK)
            {
                textBox1.Text = myOpenFileDiag.FileName;
            }

        }

        private void button1_Click(object sender, EventArgs e)
        {
            if (textBox1.Text.Trim() != String.Empty)
            {
                try
                {
                    // TODO: Validate list of tables (i.e. sheets in excel workbook)
                    string[] strTables = { "Point Survey$", "Calculation$", "Metrics$", "Diagnostics$" };

                    Form2 selectTable = new Form2(strTables);
                    selectTable.ShowDialog(this);
                    selectTable.Dispose();
                    if ((mySelectedTable != string.Empty) && (mySelectedTable != null))
                    {
                        DataTable dt = GetTable(textBox1.Text, mySelectedTable);
                        dataGridView1.DataSource = dt.DefaultView;
                    }
                }
                catch (Exception ex)
                {
                    MessageBox.Show(ex.Message.ToString());
                }
            }
            else
            {
                MessageBox.Show("Please select excel file to begin.");
            }

        }

        private DataTable GetTable(string myPath, string myTable) {
            string connectionString = "Driver={Microsoft Excel Driver (*.xls, *.xlsx, *.xlsm, *.xlsb)};DBQ=" + myPath;
            using (OdbcConnection conn = new OdbcConnection(connectionString))
            {
                conn.Open();
                string strQuery = "SELECT * FROM [" + myTable + "]";
                OdbcCommand cmd = new OdbcCommand(strQuery, conn);
                OdbcDataAdapter adapter = new OdbcDataAdapter(cmd);
                System.Data.DataSet ds = new System.Data.DataSet();
                adapter.Fill(ds);
                return ds.Tables[0];
            }
        }

    }
}
