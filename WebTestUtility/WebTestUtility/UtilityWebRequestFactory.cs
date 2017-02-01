using System;
using System.IO;
using System.Net;
using System.Security.Cryptography.X509Certificates;
using System.Threading.Tasks;

namespace WebTestUtility
{
    public interface IUtilityWebRequest
    {
        Task<string> ExecuteRequest();
    }

    public static class UtilityWebRequestFactory
    {
        public static IUtilityWebRequest CreateRequest(WebServerConfigViewModel config)
        {
            Log.Information("Factory - Creating request URI");
            Uri requestUri;
            if (!Uri.TryCreate(config.RequestUrl, UriKind.Absolute, out requestUri))
            {
                throw new Exception("No request url");
            }

            Log.Information("Factory - Setting user information");
            NetworkCredential userCred = GetUser(config);

            if (config.UseHttps)
            {
                Log.Information("Factory - set HTTPS Protocol");
                SetProtocol();

                Log.Information("Factory - Setting HTTPS certificate");
                X509Certificate2 certificate = GetClientCertificate(config);

                Log.Information("Factory - Creating HTTPS request");
                return new UtilityHttpsWebRequest(requestUri, userCred, certificate);
            }
            else
            {
                Log.Information("Factory - Creating HTTP request");
                return new UtilityHttpWebRequest(requestUri, userCred);
            }
        }

        private static void SetProtocol()
        {
            ServicePointManager.SecurityProtocol = SecurityProtocolType.Ssl3 | SecurityProtocolType.Tls | SecurityProtocolType.Tls11 | SecurityProtocolType.Tls12;
        }

        private static NetworkCredential GetUser(WebServerConfigViewModel config)
        {
            NetworkCredential userCred;
            switch (config.UserAccountOption)
            {
                case UserOption.UseSpecifiedAccount:
                    userCred = CreateSystemCredential(config);
                    break;

                case UserOption.UseCurrentAccount:
                    Log.Information("Factory - Using default credentials");
                    userCred = CredentialCache.DefaultNetworkCredentials;
                    break;

                default:
                    throw new Exception("No user credentials.");
            }

            return userCred;
        }

        private static NetworkCredential CreateSystemCredential(WebServerConfigViewModel config)
        {
            Log.Information($"Factory - Creating system credentials. User: {config.UserName}");

            if (string.IsNullOrWhiteSpace(config.UserName))
            {
                throw new Exception("No user name.");
            }
            if (config.Password.Length == 0)
            {
                throw new Exception("No password.");
            }

            return new NetworkCredential(config.UserName, config.Password);
        }

        private static X509Certificate2 GetClientCertificate(WebServerConfigViewModel config)
        {
            Log.Information($"Factory - Loading local certificate. Path: {config.CertFilePath}");

            if (string.IsNullOrWhiteSpace(config.CertFilePath))
            {
                throw new Exception("No certificate file path.");
            }

            return new X509Certificate2(File.ReadAllBytes(config.CertFilePath));
        }
    }
}